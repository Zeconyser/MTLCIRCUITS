#!/bin/bash
 
mkdir -p coregistration

sides=( "left" "right" )


for side in ${sides[@]}; do 
	output=coregistration/rASHS_${side}.nii.gz 

	if [[ -f "$output" ]]; then
        	echo "$output already exists. Skipping."
        	continue
	fi
		
	flirt -in sub-aae961_ASHS_${side}_3t_fu.nii.gz -ref t1_3t_fu.nii -applyxfm -init T2_to_T1.mat -interp nearestneighbour -out ${output} 
done	


flirt -in sub-159_ses-fu1SKYRA_eddy.nii.gz -ref t1_3t_fu.nii -out coregistration/rdwi.nii.gz -omat dwi_to_T1.mat -dof 6 -cost normmi
