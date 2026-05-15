import nibabel as nib
import numpy as np
 
 
tissue_path = r"C:\Users\z5171263\Downloads\Augusta_nifti\m2m_UA010\label_prep\tissue_labeling_upsampled.nii.gz"
tissue_img = nib.load(tissue_path)
tissue_data = tissue_img.get_fdata().astype(np.uint8)
 
 
skull_path = r"C:\Users\z5171263\Downloads\Augusta_nifti\UA010_skull_segment\Segmentation.nii"
skull_img = nib.load(skull_path)
skull_data = skull_img.get_fdata()
 
SKULL_LABEL = 7  
 
 
print("Tissue shape:", tissue_data.shape)
print("Skull shape:", skull_data.shape)
 
 
SPONGY_LABEL = 8  # spongy bone
 
tissue_data[(skull_data == 1) & (tissue_data != SPONGY_LABEL)] = SKULL_LABEL
 
 
out_img = nib.Nifti1Image(tissue_data, tissue_img.affine, tissue_img.header)
 
out_path = r"C:\Users\z5171263\Downloads\Augusta_nifti\m2m_UA010\label_prep\tissue_labeling_upsampled_UPDATED.nii.gz"
nib.save(out_img, out_path)
 
print("Saved to:", out_path)