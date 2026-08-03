
MTLCIRCUIT - PIPELINE

Abstract

Entorhinal - hippocampal circuitry contributes crucially to episodic memory function, yet is also vulnerable to aging and pathology. In this project we explore if longitudinal changes in structural connectivity between MSP and TSP assosciated subfields, are predicted by baseline AD biomarker levels (including enotrhinal Tau PET). We also investigate laminar microstructural myelination integrity in the same subfields, possible moderating effects by MTL vascularization and effects on longitudinal cognitive performance. 

Requirements 

This Pipeline uses several software suits and both 3T and 7T structural MRI.

Software (Most recent version each):
 
	FSL - for most of the registration tasks, as well as diffusion volume extraction, diffusion preprocessing (Topup and Eddy corrections) and tissue segmentation
	MRtrix - for the diffusion-based tractography pipeline and optional Diffusion Tensor 	Imaging Models. 
	ASHS - for automated segmentation of MTL subfields at 3T and 7T 
	LAYNII - for generation of Layers on MTL masks 
	(optinal: MATLAB for registrations, tissue segmenations, NORDIC denoising if 7T diffusion data becomes available) 


 
Set-up:

For the Pipeline main.sh script to work, the project folder should have the following structure: 
 
Project Folder (with subject Folders: 
	|
	 -> Subject: 
		-> 3T 
	       |    |
	       |    |-> bl / fu
	       |        - dwi_eddy.nii
	       |        - dwi.bvec	
               |        - dwi.bval 	
	       |        - T1w.nii 
	       |        - T2w.nii 
	       |        - PET	
	       |	
		-> 7T
		    |
                    |-> T1w/qT1
			T2w
			ToF		


Pipeline: 








