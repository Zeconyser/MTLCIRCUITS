import nibabel as nib
import numpy as np 
import os


base = "/mnt/d/Test_MTLCIRCUIT" 
timepoints = ["bl", "fu"]
hemis  = ["left", "right"]
pathways = ["MSP", "TSP"]
rois= {"CA1": 1, 
       "CA3" : 4,
       "DG" : 3, 
       "SUB" : 8, 
       "ERC": 9, 
       }
for t in timepoints: 
    print(f"\nProcessing {t}...")
    for hemi in hemis: 
        print(f"\nHemisphere: {hemi}...")

        ashs_path = os.path.join(base, t, "coregistration" , f"rASHS_{hemi}.nii.gz")
        ashs = nib.load(ashs_path) 
        ashs_img = ashs.get_fdata() 
        affine = ashs.affine
        seg = np.zeros_like(ashs_img)

        for roi, label in rois.items():
            print(f"\nExtracting Subfield {roi}")
            seg[np.where(ashs_img == label)] = label
            dest = os.path.join(base, t,  "subfields")
            if not os.path.exists(dest):
                print("\nCreating Destination Folder...")
                os.mkdir(dest)
            seg_file = os.path.join(dest, f"{roi}_{hemi}_{t}.nii.gz")
            print("\nSaving...")
            seg_image = nib.Nifti1Image(seg, affine) 
            nib.save(seg_image, seg_file) 




