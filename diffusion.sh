#!/bin/bash

BASE="/Volumes/LaCieJ/Test_MTLCIRCUIT/fu/Tracts"


dwi2fod -fslgrad  ${BASE}/dwi.bvec ${BASE}/dwi.bval  msmt_csd ${BASE}/dwi.mif ${BASE}/wm_response.txt ${BASE}/wm_fod.mif ${BASE}/wm_response.txt ${BASE}/gm_fod.mif ${BASE}/csf_response.txt ${BASE}/csf_fod.mif -mask ${BASE}/nodif_brain_mask.mif

