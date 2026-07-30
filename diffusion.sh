#!/bin/bash


# MRtrix3 Tractography Pipeline in native dMRI space, for Baseline and Follow Up Timepoints. 
# Expects all required images and Files to be in a single "Tracts" folder.
# Expects dwi Volumes to be fully preprocessed and bvec file to be eddy-rotated.
#
# Required Files: dwi.mif, dwi.bvecs, dwi.bvals, brain_mask.mif

# NOTE: For the Pipeline in native diffusion Space, probabilistic instead of ACT is used, as the course 2mm resolution does not seem suited to delineate a GM-WM boundary.  

# The Pipeline is then: dwi2response -> 3 tissue response function estimation
#			dwi2fod      -> Estimate 3 tissue Fiber Orientation Densities
#			tckgen	     -> Probabilistic Tractography with 10M streamlines 
#			tcksift      -> Improve tractogram fit and downsample to 2M streamlines
#			tcksift2     -> Get weighted Streamline Densities		


TIMEPOINTS=( "bl" "fu" )

for t in "${TIMEPOINTS[@]}"; do
	
	printf "\nProcessing Timepoint ${t}...\n"
	
	BASE="/Volumes/LaCieJ/Test_MTLCIRCUIT/$t/Tracts"
	DWI="${BASE}/dwi.mif"	
	
	FOD_OUT="${BASE}/wm_fod.mif"
	TCK_OUT="${BASE}/WB_PROB_10M.tck"
	
	# 1. Estimate GM, WM, CSF Response Functions
	
	if [ ! -f ${BASE}/wm_response.txt ]; then

		dwi2response \
			dhollander \
			$DWI ${BASE}/wm_response.txt ${BASE}/gm_response.txt ${BASE}/csf_response.txt \
			-fslgrad ${BASE}/dwi.bvec ${BASE}/dwi.bval

	else 
		echo "Response Functions already exist. Skipping..."
	fi
	

	# 2. Calculate GM, WM, CSF Fiber Orientation Densities

	if [ ! -f ${FOD_OUT} ]; then
		
		echo "Creating FODs..."
			  
		dwi2fod \
			 -fslgrad  ${BASE}/dwi.bvec ${BASE}/dwi.bval \
			msmt_csd ${DWI} \
		 	${BASE}/wm_response.txt ${BASE}/wm_fod.mif \
		 	${BASE}/gm_response.txt ${BASE}/gm_fod.mif \
		 	${BASE}/csf_response.txt ${BASE}/csf_fod.mif \
		 	-mask ${BASE}/brain_mask.mif

	else 
		echo "FODs already exist. Skipping..." 
	fi

	# 3. Create the probabilistic Whole-Brain Tractogram

	if [ ! -f ${TCK_OUT} ]; then

		echo "Creating WB ACT Tractogram with 10M streamlines..." 

		tckgen \
			${FOD_OUT} \
			${TCK_OUT} \
			-algorithm IFOD2 \
			-backtrack \
			-seed_dynamic ${FOD_OUT} \
			-select 10M
				
	# 4. Improve Tractogram with SIFT

		echo "Improving Tractogram with SIFT..."


		tcksift \
			${TCK_OUT} \
			${FOD_OUT} \
			"${BASE}/WB_PROB_SIFT_2M.tck" \
			-term_number 2M

		echo "Now Generating Streamline Weights with SIFT2..."

		tcksift2 \
			"${BASE}/WB_PROB_SIFT_2M.tck" \
			${FOD_OUT} \
			"${BASE}/SIFT2_WEIGHTS.txt" 
			


	else 
		echo "Tractograms already exist. Skipping..."
	fi
	
	printf "\nFinished Timepoint $t.\n"

done
