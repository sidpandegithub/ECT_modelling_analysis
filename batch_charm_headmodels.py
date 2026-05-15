import os
import subprocess
folder = r"C:\Users\z5171263\Downloads\Ramsay_Northside_nii_tryingimageJ"
for file in os.listdir(folder):
    if file.endswith(".nii"):
        subject = os.path.splitext(file)[0]
        nii = os.path.join(folder, file)
        cmd = ["charm", subject, nii]
        subprocess. run(cmd)