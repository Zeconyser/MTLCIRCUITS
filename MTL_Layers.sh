#!/bin/bash

# Base directories
RIM_BASE="/Volumes/SFB/B04/MTLCIRCUITS/Layers"
OUT_BASE="/Volumes/SFB/B04/MTLCIRCUITS/Layers"

# Path to LN2_LAYERS
LN2_LAYERS_BIN=~/laynii/LN2_LAYERS

# Number of layers
NR_LAYERS=21

# Loop over all subjects in Rim_Files
for subj_dir in "${RIM_BASE}"/sub-*; do
    subj_id=$(basename "${subj_dir}")
    echo "Processing subject: ${subj_id}"

    # Create output directory
    out_dir="${OUT_BASE}/${subj_id}/layer_files"
    mkdir -p "${out_dir}"

    # Loop over both hemispheres
    for side in left right; do
        rim_file="${subj_dir}/${subj_id}_ASHS_${side}_upsampled_rim.nii"

        if [[ -f "${rim_file}" ]]; then
            echo "  Running LN2_LAYERS for ${side}..."
            "${LN2_LAYERS_BIN}" \
                -rim "${rim_file}" \
                -nr_layers "${NR_LAYERS}" \
                -equivol
                -output "${out_dir}/${subj_id}_${side}_N21.nii"
        else
            echo "  WARNING: Rim file not found for ${side}: ${rim_file}"
        fi
    done

    echo "Finished subject: ${subj_id}"
    echo "------------------------------------"
done

echo "All subjects processed!"
