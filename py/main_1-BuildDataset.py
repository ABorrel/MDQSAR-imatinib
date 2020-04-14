import ChEMBLTable
import pathFolder



##########
#  MAIN  #
##########

# gleevec = CHEMBL941
# outlier = CHEMBL2382016

###################
# define paths directory
###################
pr_root = pathFolder.PR_ROOT
pr_result = pathFolder.PR_RESULT
pr_data = pathFolder.PR_DATA

##############
# 1. prepare dataset from ChEMBL
##############


######################
# TABLE CHEMBL  => create dataset from CHEMBL
######################
pr_tabin = pr_data + "/bioactivity-TK-ABL_CHEMBL1862.txt"
pr_result_tabclean = pathFolder.createFolder(pr_result + "CHEMBL_dataset/")
laffselected = ["IC50", "Ki"] # define affinity selected
lBAout = ["CHEMBL3705971"] # outlier
lBAout = []
ctabChEMBL = ChEMBLTable.ChEMBLTable(pr_tabin, pr_result_tabclean, laffselected, lBAout)
ctabChEMBL.CleanCHEMBLFileProtAff()# clean table
ctabChEMBL.writeTableAff() # write the affinity table only

#####
# Extract SMILES cleaned
#####
pr_chem = pathFolder.createFolder(pr_result + "SMI_chem/")
ctabChEMBL.getChemSMI(pr_chem)

##############
# Analyse dataset select
#############
ctabChEMBL.analysisTable(pr_result_tabclean)
