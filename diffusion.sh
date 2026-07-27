#!/bin/bash



dwi2fod msmt_csd dwi.mif wm_response.txt wm_fod.mif gm_response.txt gm_fod.mif csf_response.txt csf_fod.mif -mask nodif_brain_mask.nii.gz

