import pathFolder

import sys
sys.path.insert(0, "./Docking_MD/") # for window dev
import Docking_process


########################################################################################################################

##########
#  MAIN  #
##########
pr_data = pathFolder.PR_DATA
pr_result = pathFolder.PR_RESULT

##################
# 2. docking
##################
# doking is realisez with the selected chemicals using Glide software with a XP function score
# protein is prepared using Maestro on 2HYY chain A. Docking is realised on the Gleevec binding site
p_protein_prepared = pr_data + "2HYY_dock.pdb"
psdf_docking_pose = pr_data + "docking_XP-2HYY/PoseXP.sdf"
pChEMBLDataset = pr_result + "CHEMBL_dataset/tab_filtered_IC50-Ki_manualedit.csv"

###############
# analysis docking and extract poses
###############
pr_docking = pathFolder.createFolder(pr_result + "2HYY_XPdock/")

cDock = Docking_process.Docking_process(psdf_docking_pose, pr_docking)
cDock.loadSDF()
cDock.get_bestdockingscore()
cDock.plot_dockScoreVSActivity(pChEMBLDataset)
cDock.get_bestPoses()
cDock.get_TopRanking(pChEMBLDataset, 5)


