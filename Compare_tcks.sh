#!/bin/bash

# Quantitative comparison: native-space-warped-to-T1 (folder "_1") vs T1-space/ACT (folder "_2")
# tractograms, for each timepoint (bl, fu).
# Requires compare_tdi.py in the same directory (or adjust PY_SCRIPT path).

PROJECT="/Volumes/LaCieJ/Test_MTLCIRCUIT"
OUT_DIR="${PROJECT}/tck_comparisons"
TIMEPOINTS=( "bl" "fu" )
PY_SCRIPT="/Volumes/LaCieJ/MTLCIRCUITS/Tractogram_comparisons.py"   # <-- update this to wherever you saved compare_tdi.py
SUMMARY_CSV="${OUT_DIR}/tractogram_comparison_summary.csv"


mkdir -p ${OUT_DIR}

# Write summary header once
echo "timepoint,n_voxels_a,n_voxels_b,n_voxels_union,n_voxels_intersection,dice,pearson_r,pearson_p,spearman_r,spearman_p" > ${SUMMARY_CSV}

for t in "${TIMEPOINTS[@]}"; do

	printf "\nComparing tractograms for Timepoint ${t}...\n"

	# "_1" = native diffusion space, warped to T1  |  "_2" = T1-space/ACT
	NATIVE_START="${PROJECT}/${t}_1"
	T1SPACE_START="${PROJECT}/${t}_2"

	NATIVE_BASE="${NATIVE_START}/Tracts"
	T1SPACE_BASE="${T1SPACE_START}/Tracts"

	T1="${T1SPACE_START}/t1.nii"   # reference grid: the T1 used for the ACT/T1-space pipeline

	NATIVE_IN_T1_TCK="${NATIVE_BASE}/rWB_PROB_SIFT_2M.tck"
	T1SPACE_TCK="${T1SPACE_BASE}/WB_ACT_SIFT_2M.tck"

	# Save comparison outputs alongside the T1-space folder (arbitrary choice - adjust if you'd rather have a dedicated comparison folder)
	TDI_NATIVE="${OUT_DIR}/${t}_tdi_native_in_t1.nii.gz"
	TDI_T1SPACE="${OUT_DIR}/${t}_tdi_t1space_${t}.nii.gz"
	TMP_CSV="${OUT_DIR}/${t}_tdi_comparison.csv"	

	# --- 1. Streamline count / length stats for both tractograms ---
	echo "--- tckstats: native-in-T1 (${t}_1) ---"
	tckstats ${NATIVE_IN_T1_TCK}

	echo "--- tckstats: T1-space/ACT (${t}_2) ---"
	tckstats ${T1SPACE_TCK}

	# --- 2. Generate TDIs on the shared T1 grid ---
	if [ ! -f ${TDI_NATIVE} ]; then
		tckmap ${NATIVE_IN_T1_TCK} ${TDI_NATIVE} -template ${T1} -vox 1 -force
	fi

	if [ ! -f ${TDI_T1SPACE} ]; then
		tckmap ${T1SPACE_TCK} ${TDI_T1SPACE} -template ${T1} -vox 1 -force
	fi

	# --- 3. Dice + voxel-wise correlation --

	python3 ${PY_SCRIPT} ${TDI_NATIVE} ${TDI_T1SPACE} "native_in_T1" "T1space_ACT" ${TMP_CSV}

	# Append to running summary (skip header line)
	tail -n +2 ${TMP_CSV} | sed "s/^/${t},/" >> ${SUMMARY_CSV}

	printf "\nFinished comparison for Timepoint $t.\n"

done

echo -e "\nAll done. Summary written to ${SUMMARY_CSV}"
