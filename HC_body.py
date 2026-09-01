import nibabel as nib
from nibabel.processing import resample_from_to
import numpy as np
import os

p_dir_mtlrsfc = "/Volumes/SFB/B04/MTLRSFC/data/nifti/derivatives/ROIs"
p_dir_mtlcircuits = "/Volumes/SFB/B04/MTLCIRCUITS/Layers"

subs = os.listdir(p_dir_mtlcircuits)
hemis = {"LH": "left", "RH": "right"}
hc_rois = ["CA1", "CA2", "CA3", "DG"]


def create_HC_body(rois):
    body = np.zeros_like(rois[0])
    for roi in rois:
        body[roi > 0] = 1
    return body


for sub in subs:
   
    for h1, h2 in hemis.items():
    
        p_out = os.path.join(
                p_dir_mtlcircuits,
                sub,
                f"{sub}_ASHS_HB_{h2}_upsampled.nii.gz"
                )
        if os.path.exists(p_out):
            print(f"skipping {sub}...Already exists.")
            continue
        
                
        print(f"\nProcessing {sub}'s {h2} side...")

        rois = []

        for roi in hc_rois:
            p_roi = os.path.join(
                p_dir_mtlrsfc,
                sub,
                "func_ashs_ROIs",
                f"r{sub}_{roi}_{h1}_HB_mask.nii"
            )

            i_roi = nib.load(p_roi)
            d_roi = i_roi.get_fdata()
            aff_roi = i_roi.affine

            rois.append(d_roi)

        hc_body = create_HC_body(rois)

        p_hi_ashs = os.path.join(
            p_dir_mtlcircuits,
            sub,
            f"{sub}_ASHS_{h2}_upsampled.nii.gz"
        )

        i_hi_ashs = nib.load(p_hi_ashs)
        d_hi_ashs = i_hi_ashs.get_fdata()
        aff_hi_ashs = i_hi_ashs.affine

        hi_body_img = nib.Nifti1Image(
            hc_body,
            aff_roi
        )

        hi_body_ashs = resample_from_to(
            hi_body_img,
            i_hi_ashs,
            order=0
        ).get_fdata()

        body_ashs = d_hi_ashs.copy()

        labels = [1, 2, 3, 4]
        hc_mask = np.isin(d_hi_ashs, labels)

        body_ashs[hc_mask & (hi_body_ashs == 0)] = 0

        i_out = nib.Nifti1Image(
            body_ashs,
            aff_hi_ashs
        )

        nib.save(i_out, p_out)

    print("Finished...Next")
