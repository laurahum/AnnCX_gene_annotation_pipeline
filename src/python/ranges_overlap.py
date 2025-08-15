#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 01:23:25 2024

@author: lahumada
Modified code from https://stackoverflow.com/questions/32480423/how-to-check-if-a-range-is-a-part-of-another-range-in-python-3-x) (Accessed in May 2022)
"""


# This script finds out whether a gene range is within another gene range
# which is necessary to filter the annotation results from BLAST, GMAP and
# Exonerate 

# Main function: find_ranges_overlap(range1, range2)

# Function to find if range 1 overlaps with range 2
def find_ranges_overlap(range1, range2):
    """
    Determine if range1 overlaps with range2.
    Overlap is defined as portion of range1 intersecting with range2.
    This function is necessary to filter the results from BLAST, GMAP, GeneWise and Exonerate
    If range 1 overlaps with range 2 = True
    """
    # Control statements
    if not range1:
        return True  # empty range is a subset of anything
    if not range2:
        return False  # non-empty range can't be subset of empty range
    if len(range1) > 1 and range1.step % range2.step:
        return False  # must have a single value or integer multiple step

    start1, end1 = range1.start, range1[-1]
    start2, end2 = range2.start, range2[-1]
    return start1 < end2 and end1 > start2
