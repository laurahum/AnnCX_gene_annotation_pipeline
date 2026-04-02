#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 19:43:13 2026

@author: lahumada
"""

"""
Convert miniprot GFF3 output to extract PAF score info from miniprot GFF3 and add to mRNA attributes.

(1) Parses ##PAF lines to extract alignment scores/tags
(2) Matches PAF to corresponding mRNA by Target/ID and coordinates
(3) Adds PAF fields to mRNA attributes
- AS:i     = Alignment score (higher = better)
- ms:i     = Miniprot internal matching score  
- np:i     = Number of sequence pieces (more = more fragments)
- fs:i     = Number of full sequence matches
- st:i     = Number of suboptimal hits (more = more alternatives)
- da:i     = Deletion count
- do:i     = Insertion count  
- mapq:i   = Mapping quality (0-255; 255=missing)
- nmatch:i = Number of matching bases
- alen:i   = Alignment block length (matches + gaps)
- qlen:i   = Query (protein) length
- cigar    = CIGAR string (alignment blocks)
- cs       = Compression string (miniprot-specific)
(4) Writes reformatted GFF3 with complete gene models + PAF info in mRNA

Main function: format_miniprot_output_to_gff_all_files(input_dir, output_dir)
"""


import os
import sys
import pandas as pd
import re


def parse_paf_comment_line(line: str):
    """
    Parse a '##PAF' line into a dict of core fields + tags.
    Expects: ##PAF <12 core fields> <tag:...:...> ...
    Returns None if the line is not a valid PAF line.
    """
    line = line.rstrip("\n")
    if not line.startswith("##PAF"):
        return None

    # Drop leading '##PAF'
    fields = line.split()[1:]
    
    core = fields[:12]
    tags = fields[12:]

    paf = {
        "qname": core[0],
        "qlen": int(core[1]),
        "qstart": int(core[2]),
        "qend": int(core[3]),
        "strand": core[4],
        "tname": core[5],
        "tlen": int(core[6]),
        "tstart": int(core[7]),
        "tend": int(core[8]),
        "nmatch": int(core[9]),
        "alen": int(core[10]),
        "mapq": int(core[11]),
    }

    # Optional tags (AS:i:..., cs:Z:..., cg:Z:..., etc.)
    for tag in tags:
        if ":" not in tag:
            continue
        # Split into 3 parts: tag, type, value
        parts = tag.split(":", 2)
        if len(parts) == 3:
            key, _typ, value = parts
            paf[key] = value

    return paf


def extract_target_name_from_attributes(attributes: str):
    """
    Extract 'Target=' value from an attributes string.
    Uses a conservative pattern (no spaces), customize if needed.
    """
    m = re.search(r"Target=([^;\s]+)", attributes)
    return m.group(1)


def read_gff3_to_dataframe(filepath: str):
    """
    Read GFF3 into:
      - header_lines: list of lines starting with '#'
      - df: pandas DataFrame with columns 0..8 (for GFF rows)
      - raw_lines: list of all lines (including headers) for index mapping
    We preserve line order by keeping row_index (position in file).
    """
    header_lines = []
    data_rows = []
    raw_lines = []

    with open(filepath, "r") as f:
        for idx, line in enumerate(f):
            raw_lines.append(line.rstrip("\n"))
            if line.startswith("#"):
                header_lines.append(line.rstrip("\n"))
                # Keep header (including ##PAF) in the DataFrame
                # but mark their type specially.
                cols = [None] * 9
                cols[0] = "HEADER"   # pseudo seqid label
                cols[2] = "HEADER"   # pseudo type
                cols[8] = line.rstrip("\n")  # store full header line
            else:
                # Normal GFF row, split into 9 columns
                parts = line.rstrip("\n").split("\t")
                cols = parts[:9]

            data_rows.append([idx] + cols)

    df = pd.DataFrame(
        data_rows,
        columns=[
            "row_index",  # original line index in file
            "seqid",
            "source",
            "type",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
        ],
    )

    return header_lines, df, raw_lines


def write_gff3_from_dataframe(df: pd.DataFrame, output_path: str):
    """
    Write a DataFrame with columns like above back to GFF3.
    HEADER rows are written as their stored attributes column
    (contains the full header line).
    **EXCLUDES ##PAF lines from output**
    """
    lines_out = []
    df_sorted = df.sort_values("row_index")

    for _, row in df_sorted.iterrows():
        # **SKIP ##PAF lines** - check attributes content directly
        if (row["seqid"] == "HEADER" and 
            row["type"] == "HEADER" and 
            str(row["attributes"]).startswith("##PAF")):
            continue  # Skip PAF lines
            
        if row["seqid"] == "HEADER" and row["type"] == "HEADER":
            # Write other header lines (##gff-version 3, ##FASTA, etc.)
            lines_out.append(str(row["attributes"]))
        else:
            cols = [
                row["seqid"],
                row["source"],
                row["type"],
                row["start"],
                row["end"],
                row["score"],
                row["strand"],
                row["phase"],
                row["attributes"],
            ]
            # Ensure everything is string and join with tabs
            cols = ["" if pd.isna(c) else str(c) for c in cols]
            lines_out.append("\t".join(cols))

    with open(output_path, "w") as f:
        f.write("\n".join(lines_out) + "\n")


def format_miniprot_gff_with_paf_scores(input_path: str, output_path: str):
    """
    New logic:
      1. Read GFF3 into a DataFrame.
      2. Find indices of '##PAF' lines.
      3. For each PAF block (from PAF row to row before next PAF / EOF):
         - Parse PAF info from that header line.
         - Find mRNA rows in this block.
         - Add PAF scores to mRNA attributes.
      4. Write the full DataFrame back to GFF3.
    """
    header_lines, df, raw_lines = read_gff3_to_dataframe(input_path)

    # Identify rows that are exactly '##PAF...' in the original text.
    # Look at HEADER rows whose 'attributes' starts with '##PAF'.
    paf_mask = (df["seqid"] == "HEADER") & df["attributes"].astype(str).str.startswith("##PAF")
    paf_indices = df.index[paf_mask].tolist()

    # Slice with row_index in the file to keep order consistent
    # But ue df.index positions for block boundaries.
    # Each block: [paf_idx, next_paf_idx) in df.index space.
    paf_indices_sorted = sorted(paf_indices)

    for i, paf_df_idx in enumerate(paf_indices_sorted):
        # Block start: current PAF row
        block_start_idx = paf_df_idx

        # Block end: row before next PAF row index, or end of df
        if i + 1 < len(paf_indices_sorted):
            next_paf_df_idx = paf_indices_sorted[i + 1]
            block_end_idx = next_paf_df_idx - 1
        else:
            block_end_idx = df.index[-1]

        # Sub-DataFrame view for this block
        block_df = df.loc[block_start_idx:block_end_idx]

        # Parse PAF from the PAF row text
        paf_line = block_df.loc[paf_df_idx, "attributes"]
        paf_info = parse_paf_comment_line(paf_line)

        # Build PAF attribute string pieces to insert into mRNA
        paf_attrs = [
            f"AS={paf_info.get('AS', '.')}",
            f"ms={paf_info.get('ms', '.')}",
            f"np={paf_info.get('np', '.')}",
            f"fs={paf_info.get('fs', '.')}",
            f"st={paf_info.get('st', '.')}",
            f"da={paf_info.get('da', '.')}",
            f"do={paf_info.get('do', '.')}",
            f"mapq={paf_info['mapq']}",
            f"nmatch={paf_info['nmatch']}",
            f"alen={paf_info['alen']}",
            f"qlen={paf_info['qlen']}",
        ]
        if "cg" in paf_info:
            paf_attrs.append(f"cigar={paf_info['cg']}")
        if "cs" in paf_info:
            paf_attrs.append(f"cs={paf_info['cs']}")

        paf_attrs_str = ";".join(paf_attrs)

        # Find mRNA rows in this block and update attributes
        # If multiple mRNA rows exist, apply to all or adjust logic as needed.
        mRNA_mask = block_df["type"] == "mRNA"
        mRNA_indices = block_df.index[mRNA_mask].tolist()
        for mrna_idx in mRNA_indices:
            old_attr = str(df.at[mrna_idx, "attributes"])

            # Append PAF attributes; ensure we keep any existing attributes first
            if old_attr.endswith(";") or old_attr == "":
                new_attr = f"{old_attr}{paf_attrs_str}"
            else:
                new_attr = f"{old_attr};{paf_attrs_str}"
            df.at[mrna_idx, "attributes"] = new_attr

    # Finally, write the modified DataFrame back to GFF3
    write_gff3_from_dataframe(df, output_path)


def format_miniprot_output_to_gff_all_files(input_dir: str, output_dir: str):
    """
    Apply the DataFrame-based PAF → GFF3 formatting to all files in a directory.
    """
    os.makedirs(output_dir, exist_ok=True)

    for filename in os.listdir(input_dir):
        output_name = filename.rsplit('.', 1)[0] + '_FORMATTED.gff3'
        Filepath_input = os.path.join(input_dir, filename)
        Filepath_output = os.path.join(output_dir, output_name)
            
        format_miniprot_gff_with_paf_scores(Filepath_input, Filepath_output)
        print(f"Saved: {output_name}")
