#!/bin/bash

# Define paths and variables
REPO_DIR="$(pwd)"
ENV_YML_FILE="${REPO_DIR}/envs/environment.yml"
PYTHON_MAIN_SCRIPT="${REPO_DIR}/AnnCX.py"
PYTHON_PRED2REF_SCRIPT="${REPO_DIR}/src/identify_pred2ref.py"
PYTHON_REARRANGEMENT_SCRIPT="${REPO_DIR}/src/identify_rearrangements.py"
PYTHON_ANNOTATION2FASTA_SCRIPT="${REPO_DIR}/src/annotation2fasta.py"

# Step 1: Create Conda environment from YAML file
echo "Creating Conda environment from ${ENV_YML_FILE}..."
conda env create -f "${ENV_YML_FILE}" || {
    echo "Error: Failed to create Conda environment from ${ENV_YML_FILE}."
    exit 1
}
echo "Conda environment created successfully."

# Step 2: Give execution rights
# Main Python script
echo "Grant execution rights to ${PYTHON_MAIN_SCRIPT}..."
chmod +x "${PYTHON_MAIN_SCRIPT}" || {
    echo "Error: Failed to grant execution rights to ${PYTHON_MAIN_SCRIPT}."
    exit 1
}
echo "Grant execution rights to ${PYTHON_PRED2REF_SCRIPT}..."
chmod +x "${PYTHON_PRED2REF_SCRIPT}" || {
    echo "Error: Failed to grant execution rights to ${PYTHON_PRED2REF_SCRIPT}."
    exit 1
}
echo "Grant execution rights to ${PYTHON_REARRANGEMENT_SCRIPT}..."
chmod +x "${PYTHON_REARRANGEMENT_SCRIPT}" || {
    echo "Error: Failed to grant execution rights to ${PYTHON_REARRANGEMENT_SCRIPT}."
    exit 1
}
echo "Grant execution rights to ${PYTHON_ANNOTATION2FASTA_SCRIPT}..."
chmod +x "${PYTHON_ANNOTATION2FASTA_SCRIPT}" || {
    echo "Error: Failed to grant execution rights to ${PYTHON_ANNOTATION2FASTA_SCRIPT}."
    exit 1
}
echo "Execution rights granted successfully."


# Final message
echo "Pipeline installation completed successfully."

