#!/bin/bash


# MRtrix3 Tractography Pipeline in native dMRI space, for Baseline and Follow Up Timepoints. 
# Expects all required images and Files to be in a single "Tracts" folder.
# Expects dwi Volumes to be fully preprocessed and bvec file to be eddy-rotated.
#
# Required Files: dwi.nii, dwi.mif, dwi.bvecs, dwi.bvals, brain_mask.mif, T1.nii, mean_b0.nii

# NOTE: For the Pipeline in T1 anatomical  Space, ACT is used as the finer 0.8mm resolution is more suited to delineate a GM-WM boundary.  

# The Pipeline is then: 
#			5ttgen       -> 5 Tissue segmentation on T1
#			5tt2gmwmi    -> Create GM-WM-seed used for ACT
# 			flirt 	     -> Register dwi + brain_mask to T1 space	
#			dwi2response -> 3 tissue response function estimation
#			dwi2fod      -> Estimate 3 tissue Fiber Orientation Densities
#			tckgen	     -> Probabilistic Tractography with 10M streamlines 
#			tcksift      -> Improve tractogram fit and downsample to 2M streamlines
#			tcksift2     -> Get weighted Streamline Densities		


TIMEPOINTS=( "bl_2" "fu_2" )

for t in "${TIMEPOINTS[@]}"; do
	
	printf "\nProcessing Timepoint ${t}...\n"
	
	START="/Volumes/LaCieJ/Test_MTLCIRCUIT/$t"
	BASE="/Volumes/LaCieJ/Test_MTLCIRCUIT/$t/Tracts"
	DWI="${START}/dwi.nii.gz"	
	T1="${START}/t1.nii"
	B0="${START}/b0_mean.nii.gz"
	
	

	# 1 Preparations: 
	  # 1.1 Create 5TT mask and gmwmseed 
	
 	if [ ! -f ${BASE}/GMWMseed.mif ]; then
		
		printf "\nCreating GM-WM Seed\n"
		
		5ttgen \
			fsl \
			$T1 \
			${BASE}/5TT.mif				          
		
		5tt2gmwmi \
			${BASE}/5TT.mif \
			${BASE}/GMWMseed.mif

	else 
                echo "\nGMWMseed already exists. Skipping..."
        fi


	  # 1.2 Register dwi to T1 and rotate bvecs
	
	if [ ! -f ${BASE}/rdwi.mif ]; then
	
		printf "\nCoregistering DWI and brain mask to T1 space..."
		
		echo "estimating transform..." 

		flirt \
			-in $B0 \
			-ref $T1 \
			-omat ${BASE}/dwi_to_t1.mat \
			-dof 6 \
			-cost normmi 
		
		echo "Rotating b-vectors..." 

		fdt_rotate_bvecs \
			${BASE}/dwi.bvec \
            ${BASE}/rdwi.bvec \
            ${BASE}/dwi_to_t1.mat
			      	
	 	
		echo "Applying Registration to dwi volume and brain mask..."
		
		flirt \
			-in $DWI \
			-ref $T1 \
			-applyxfm \
			-init ${BASE}/dwi_to_t1.mat \
			-out ${BASE}/rdwi.nii.gz

		flirt \
			-in ${BASE}/brain_mask.nii.gz \
			-ref $T1 \
			-applyxfm \
			-init ${BASE}/dwi_to_t1.mat \
            -interp nearestneighbour \
			-out ${BASE}/rbrain_mask.nii.gz



		echo "Converting to .mif..." 

		

		
		mrconvert ${BASE}/rdwi.nii.gz ${BASE}/rdwi.mif
		mrconvert ${BASE}/rbrain_mask.nii.gz ${BASE}/rbrain_mask.mif

	else 
                echo "\nDWI and b-vectors already in T1 space..."
        fi


	
	FOD_OUT="${BASE}/wm_fod.mif"
	TCK_OUT="${BASE}/WB_PROB_10M.tck"
	DWI="${BASE}/rdwi.mif" 
	
	# 2. Estimate GM, WM, CSF Response Functions
	
	if [ ! -f ${BASE}/wm_response.txt ]; then

		dwi2response \
			dhollander \
			$DWI ${BASE}/wm_response.txt ${BASE}/gm_response.txt ${BASE}/csf_response.txt \
			-fslgrad ${BASE}/rdwi.bvec ${BASE}/dwi.bval

	else 
		echo "Response Functions already exist. Skipping..."
	fi
	

	# 3. Calculate GM, WM, CSF Fiber Orientation Densities

	if [ ! -f ${FOD_OUT} ]; then
		
		echo "Creating FODs..."
			  
		dwi2fod \
			 -fslgrad  ${BASE}/rdwi.bvec ${BASE}/dwi.bval \
			msmt_csd ${DWI} \
		 	${BASE}/wm_response.txt ${BASE}/wm_fod.mif \
		 	${BASE}/gm_response.txt ${BASE}/gm_fod.mif \
		 	${BASE}/csf_response.txt ${BASE}/csf_fod.mif \
		 	-mask ${BASE}/rbrain_mask.mif

	else 
		echo "FODs already exist. Skipping..." 
	fi

	# 4. Create the probabilistic Whole-Brain Tractogram

	if [ ! -f ${TCK_OUT} ]; then

		echo "Creating WB ACT Tractogram with 10M streamlines..." 

		tckgen \
			${FOD_OUT} \
			${TCK_OUT} \
			- act ${BASE}5TT.mif \
			-backtrack \
			-seed_gmwmi ${BASE}/GMWMseed.mif \
			-select 10M
				
	# 5. Improve Tractogram with SIFT

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
			"${BASE}/SIFT2_WEIGHTS.txt" \
			-act 5TT.mif
			


	else 
		echo "Tractograms already exist. Skipping..."
	fi
	
	printf "\nFinished Timepoint $t.\n"

done
