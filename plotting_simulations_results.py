# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 09:41:35 2026

@author: z5171263
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TP_sim =  pd.read_excel("C:/Users/z5171263/Downloads/Northside_TP_sims/assemblynet_TP_sim_ROI_summary.xls")
FP_sim =  pd.read_excel("C:/Users/z5171263/Downloads/Northside_FP_sims/assemblynet_FP_sim_ROI_summary.xls")

TP_sim = TP_sim.drop(columns=['Row'])
FP_sim = FP_sim.drop(columns=['Row'])

##########
mean_FP = FP_sim.mean()
std_FP  = FP_sim.std()

mean_TP = TP_sim.mean()
std_TP  = TP_sim.std()


rois = mean_FP.index
x = np.arange(len(rois))
width = 0.35

plt.figure()

plt.bar(x - width/2, mean_FP.values, width, yerr=std_FP.values, label='FP')
plt.bar(x + width/2, mean_TP.values, width, yerr=std_TP.values, label='TP')

plt.xticks(x, rois, rotation=90)
plt.ylabel("Mean E-field")
plt.title("ROI Comparison (FP vs TP)")
plt.legend()

plt.tight_layout()
plt.show()

###############
# Build base ROI names (LH_ → strip prefix)
base_rois = [col[3:] for col in mean_FP.index if col.startswith('LH_')]

# Means
LH_FP = mean_FP[['LH_' + roi for roi in base_rois]]
RH_FP = mean_FP[['RH_' + roi for roi in base_rois]]

LH_TP = mean_TP[['LH_' + roi for roi in base_rois]]
RH_TP = mean_TP[['RH_' + roi for roi in base_rois]]

# STDs
LH_FP_std = std_FP[['LH_' + roi for roi in base_rois]]
RH_FP_std = std_FP[['RH_' + roi for roi in base_rois]]

LH_TP_std = std_TP[['LH_' + roi for roi in base_rois]]
RH_TP_std = std_TP[['RH_' + roi for roi in base_rois]]

# Correct x for THIS plot
x = np.arange(len(base_rois))
width = 0.2

plt.figure(figsize=(20, 8))

# FP
plt.bar(x - 1.5*width, LH_FP.values, width,
        yerr=LH_FP_std.values, capsize=4, label='LH (FP)')

plt.bar(x - 0.5*width, RH_FP.values, width,
        yerr=RH_FP_std.values, capsize=4, label='RH (FP)')

# TP
plt.bar(x + 0.5*width, LH_TP.values, width,
        yerr=LH_TP_std.values, capsize=4, label='LH (TP)')

plt.bar(x + 1.5*width, RH_TP.values, width,
        yerr=RH_TP_std.values, capsize=4, label='RH (TP)')

plt.xticks(x, base_rois, rotation=90)
plt.ylabel("Mean E-field")
plt.title("ROI Comparison (FP vs TP, LH vs RH)")
plt.legend()

plt.tight_layout()
plt.savefig("figure.png")
plt.show()




####################################################

# Function to reshape + summarise
def summarise(df, condition):
    # Convert wide → long
    long_df = df.melt(var_name='ROI_full', value_name='Efield')
    
    # Extract hemisphere + base ROI
    long_df['Hemisphere'] = long_df['ROI_full'].str[:2]
    long_df['ROI'] = long_df['ROI_full'].str[3:]
    
    long_df['Condition'] = condition
    
    # Group + compute stats
    summary = (long_df
               .groupby(['Condition', 'Hemisphere', 'ROI'])
               .agg(
                   Mean=('Efield', 'mean'),
                   Std=('Efield', 'std'),
                   N=('Efield', 'count')
               )
               .reset_index())
    
    summary['SEM'] = summary['Std'] / np.sqrt(summary['N'])
    
    return summary

summary_FP = summarise(FP_sim, 'FP')
summary_TP = summarise(TP_sim, 'TP')

summary_df = pd.concat([summary_FP, summary_TP], ignore_index=True)


summary_df.to_excel("C:/Users/z5171263/OneDrive - UNSW/Desktop/analysis/master_sheets/Volbrain_sheets_Northside18only/Northside_VolBrain_summary.xlsx", index=False)


