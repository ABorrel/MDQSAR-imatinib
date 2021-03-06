import pathFolder

import sys
sys.path.insert(0, "./Docking_MD/")
import MD_run
import MD_globalAnalysis


##########
#  MAIN  #
##########
pr_data = pathFolder.PR_DATA
pr_result = pathFolder.PR_RESULT

##################
# 3. Molecular Dynamics
##################
# MD are realised using Schrodinger GDesmond.
# for each MD ligand by frame, binding sites and proteins
# and compute RMSD after a superimposition using TMalign
p_protein_prepared = pr_data + "2HYY_MD.pdb" # the protein have to be prepared before
pr_pose_selected = pr_result + "2HYY_XPdock/BEST_POSES/"
pChEMBLDataset = pr_result + "CHEMBL_dataset/tab_filtered_IC50-Ki_manualedit.csv"

# folder with all of the MD
pr_MDrun = pathFolder.createFolder(pr_result + "MD-run/")

# folder with frames for next step
pr_MDout = pathFolder.createFolder(pr_result + "MD-out/")

#########
# Parameter MD
#########
timeMD = "15000.0"
timeframe = "10.0"
stepWait = 9
nbGPU = 3 
nbCPU = 10
stepFrame = 1 # reduce the number of extracted frames
nameLig = "UNK"
water = 0
BSCutoff = 6.0 # Binding site 


#########
# RUN MD for the set of poses and compute RMSD
#########
#MD_run.computeMD(pr_pose_selected, pr_MDrun, p_protein_prepared, pr_MDout, nameLig, BSCutoff, timeMD, timeframe, stepWait, stepFrame, water, nbCPU, nbGPU)


###########
# analysis MD quality
###########
pr_globalQuality = pathFolder.createFolder(pr_result + "MD-qualityCheck/")
cMDglobalQuality = MD_globalAnalysis.MD_globalAnalysis(pChEMBLDataset, pr_MDrun, pr_MDout, pr_globalQuality)
cMDglobalQuality.buildRMDSSheets()
cMDglobalQuality.mergeRMSDSheets()
# next merge in one pdf for publication




