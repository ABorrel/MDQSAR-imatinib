import ChEMBLTable
#import ligand
import pathFolder
#import runExternalSoft
#import MCS
#import parseSDF
#import MD
#import FPI
#import PDB
#import cpdClustering
#import MDdescriptors
#import QSARModeling
#import analyzeDBLig
#import builderDescMatrix
#import clusterDesc
#import dockingAnalysis
#import MDanalysis

#from os import listdir, makedirs
#from re import search
#from time import sleep


#########################
##   MAIN FUNCTIONS    ##
#########################

#def FPIMatrix(cdocking, pprotein, prFPI):

#    pmatrixFPI = prFPI + "MFPI.txt"

#    i = 0
#    imax = len(cdocking.lposefiles)
    #imax = 3
#    while i < imax:
#        j = i + 1
#        while j < imax:
#            cprot = PDB.PDB(PDB_input=pprotein, hydrogen=1)
#            pligPDBi = runExternalSoft.babelConvertSDFtoPDB(cdocking.lposefiles[i])
#            pligPDBj = runExternalSoft.babelConvertSDFtoPDB(cdocking.lposefiles[j])

#            print pligPDBi
#            cposei = PDB.PDB(PDB_input=pligPDBi, hydrogen=1)
#            cposej = PDB.PDB(PDB_input=pligPDBj, hydrogen=1)

#            sFPIi = FPI.ligFPI(cPDB=cprot, ligin=cposei, prFPI=prFPI)
#            sFPIj = FPI.ligFPI(cPDB=cprot, ligin=cposej, prFPI=prFPI)

#            dout = sFPIi.compareFPI(sFPIj)

#            print dout

#            j = j + 1
#        i = i + 1


#def dockingAnalysis(psdfDoking, ltableCpd, ptableCpd, prpose, pranalysis):

#    sdocking = parseSDF.sdf(psdfDoking)
#    sdocking.parseSDF()
#    sdocking.splitPoses(prpose)

#    dscore = sdocking.get_dockingscore()
#    dockingScoreAnalysis(dscore, ltableCpd, ptableCpd, pranalysis)
#    return dscore



#def computeMD(prLigand, prMD, pprotein, pranalysis, nameLig, BSCutoff, timeMD, timeframe, stepWait, stepFrame, water, nbCPU, nbGPU):
#    """
#    :param prLigand: folder of ligands or poses
#    :param pprotein: protein in PDB already prepared using the maestro
#    :param timeMD: time MD in ps
#    :param timeframe: time between two frame
#    :param stepWait: time to wait for multiprocess
#    :param stepFrame: frame extraction
#    :param water: extract water molecules
#    :param nbCPU: CPU max in parrallel
#    :param nbGPU: GPU max in paralele, random selection
#    :return:
#    """
#    # 1. Merge poses and proteins
#    cMDs = MD.MD(prMD, pranalysis, water, timeMD, timeframe, stepWait, nbGPU, nbCPU, stepFrame)
#    cMDs.initialisation(prLigand, pprotein)
#    cMDs.runMultipleMD()# run MD

    # 2. Extract frames
    # extract frame
#    cMDs.centerFrame()
#    cMDs.extractFrame()

    # extract BS and ligand
#    cMDs.extractLigBSbyFrame(BSCutoff, nameLig, clean=0)

    # for complete run have to control all job done
    #sleep(100)
    # 3. compute RMSD
#    cMDs.analyseRMSD()



#def computeMDdesc(prMD, prout, istart=0, iend=0, descLig=1, descBS=1, descFPI=1 ):

#    lpMDresult = listdir(prMD)
#    if iend == 0:
#        iend = len(lpMDresult)
#    print prMD

#    i = istart
#    for MDresult in lpMDresult[istart:iend]:
#        jobname = MDresult.split("_")[0]
#        prlig = prMD + MDresult + "/lig/"
#        prBS = prMD + MDresult + "/BSs/"
#        prframes = prMD + MDresult + "/framesMD/"
#        prMDdesc = pathFolder.createFolder(prout + jobname + "/")


        # control run
#        print jobname, i
#        i +=1

        # compute different descriptors
#        if not "cMD" in locals().keys():
#            cMD = MDdescriptors.MDdescriptors(jobname, prlig, prBS, prframes, prMDdesc)
#        else:
#            cMD.jobname = jobname
#            cMD.prlig = prlig
#            cMD.prBSs = prBS
#            cMD.prframe = prframes
#            cMD.prout = prMDdesc

#        if descLig == 1:
#            cMD.computeLigDesc()
#        if descBS == 1:
#            cMD.computeBSDesc()
#        if descFPI == 1:
#            cMD.computeFPI()




########################################################################################################################
########################################################################################################################

##########
#  MAIN  #
##########

# gleevec = CHEMBL941
# outlier = CHEMBL2382016

###################
# define path  
###################
pr_root = "./../../"
pr_result = pathFolder.createFolder(pr_root + "results/")
pr_data = pathFolder.createFolder(pr_root + "data/")

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
#ctabAll.analysisTable(prCHEMBL)


################
# Extract chemical
################
pr_chem = pathFolder.createFolder(pr_result + "SMI_chem/")
ctabAll.extractChemSMI(pr_chem)
ddd





prHome = pathFolder.PR_HOME
prCHEMBL = pathFolder.createFolder(pathFolder.PR_RESULT + "CHEMBL/")
pCHEMBL = prHome + "imatinib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862.txt"
pCHEMBLout = prHome + "imatinib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862_filtered.txt"
laffselected = ["IC50", "Ki", "Kd"]
lBAout = ["CHEMBL3705971"]
lBAout = []


ctabAll = ChEMBLTable.CleanCHEMBLFileProtAff(pCHEMBL, pCHEMBLout, laffselected, lBAout)
paff = ctabAll.writeTableAff(prCHEMBL + "AffAllcurated")
#ctabAll.analysisTable(prCHEMBL)

sss

#####################################
#####################################
#  Docking analysis -SP-XP docking  #
#####################################
#####################################


# SP-2HYY #
###########
#psdfDokingSP = "/home/aborrel/imitanib/docking/dockingSP_2hyy/dockingpose.sdf"
#prDockingPoseSP = "/home/aborrel/imitanib/results/dockingposeSP_2HYY/"
#pprotein_2HYY = "/home/aborrel/imitanib/protein/2HYY_dock.pdb"
#pranalysis_SP_2HYY = pathFolder.analyses("2HYY_SPdock")

# analysis
#dockingAnalysis(psdfDokingSP, ltab, pCHEMBLClean, prDockingPoseSP, pranalysis_SP_2HYY)

# XP #
######
# native - 2HYY #
psdfDokingXP_2HYY = prHome + "imatinib/docking/dockingXP_2hyy/PoseXP.sdf"
prDockingPoseXP_2HYY = prHome + "imatinib/results/dockingbestposeXP_2HYY/"
pprotein_2HYY = prHome + "imatinib/protein/2HYY_dock.pdb"
pranalysis_XP_2HYY = pathFolder.analyses("2HYY_XPdock")
prMDanalysis = prHome + "imatinib/results/MDanalysis/"

# analysis
# 1. poses
#sdocking = parseSDF.sdf(psdfDokingXP_2HYY, prDockingPoseXP_2HYY)
#sdocking.parseSDF()
#sdocking.get_dockingscore()
#sdocking.get_bestPose()# based on glide score

#2. Scores correlation scores
#dockingAnalysis.dockingScoreAnalysis(sdocking.docking, ctabAll.table, pCHEMBLout, pranalysis_XP_2HYY)

#3. Top chemical
#nbrank = 900
#prrank = pathFolder.createFolder(pranalysis_XP_2HYY + "Rank/")
#dockingAnalysis.rankingTop(sdocking, ctabAll.table, prDockingPoseXP_2HYY, prrank, nbrank)

#4. correlation plot with RMSD
#prRMSD = pathFolder.createFolder(pranalysis_XP_2HYY + "MD/")
#dockingAnalysis.plotRMSDVSDockingScore(sdocking.docking, ctabAll.table, pCHEMBLout, prMDanalysis, prRMSD)


###############################
################################
####### MUTATED PROTEINS  ######
################################

# mutated 3QRJ
#psdfDokingXP_3QRJ = "/home/aborrel/imitanib/docking/dockingXP_3QRJ/3QRJ_XP_pose.sdf"
#prDockingPoseXP_3QRJ = "/home/aborrel/imitanib/results/dockingposeXP_3QRJ/"
#pprotein_3QRJ = "/home/aborrel/imitanib/protein/3QRJ_dock.pdb"
#pranalysis_XP_3QRJ = pathFolder.analyses("3QRJ_XPdock")

# analysis
#dockingAnalysis(psdfDokingXP_3QRJ, ltab, pCHEMBLClean, prDockingPoseXP_3QRJ, pranalysis_XP_3QRJ)


# mutated 2F4J
#psdfDokingXP_2F4J = "/home/aborrel/imitanib/docking/dockingXP_2F4J/2F4J_XP_pose.sdf"
#prDockingPoseXP_2F4J = "/home/aborrel/imitanib/results/dockingposeXP_2F4J/"
#pprotein_2F4J = "/home/aborrel/imitanib/protein/2F4J_dock.pdb"
#pranalysis_XP_2F4J = pathFolder.analyses("2F4J_XPdock")

# analysis
#dockingAnalysis(psdfDokingXP_2F4J, ltab, pCHEMBLClean, prDockingPoseXP_2F4J, pranalysis_XP_2F4J)


# mutated 2FO0
#psdfDokingXP_2FO0 = "/home/aborrel/imitanib/docking/dockingXP_2FO0/2FO0_XP_pose.sdf"
#prDockingPoseXP_2FO0 = "/home/aborrel/imitanib/results/dockingposeXP_2FO0/"
#pprotein_2FO0 = "/home/aborrel/imitanib/protein/2FO0_dock.pdb"
#pranalysis_XP_2FO0 = pathFolder.analyses("2FO0_XPdock")

# analysis
#dockingAnalysis(psdfDokingXP_2FO0, ltab, pCHEMBLClean, prDockingPoseXP_2FO0, pranalysis_XP_2FO0)

# mutated 2G2H
#psdfDokingXP_2G2H = "/home/aborrel/imitanib/docking/dockingXP_2G2H/2G2H_XP_pose.sdf"
#prDockingPoseXP_2G2H = "/home/aborrel/imitanib/results/dockingposeXP_2G2H/"
#pprotein_2G2H = "/home/aborrel/imitanib/protein/2G2H_dock.pdb"
#pranalysis_XP_2G2H = pathFolder.analyses("2G2H_XPdock")

# analysis
#dockingAnalysis(psdfDokingXP_2G2H, ltab, pCHEMBLClean, prDockingPoseXP_2G2H, pranalysis_XP_2G2H)


# mutated 2GQG - verif
#psdfDokingXP_2GQG = "/home/aborrel/imitanib/docking/dockingXP_2GQG/2GQG_XP_pose.sdf"
#prDockingPoseXP_2GQG = "/home/aborrel/imitanib/results/dockingposeXP_2GQG/"
#pprotein_2GQG = "/home/aborrel/imitanib/protein/2GQG_dock.pdb"
#pranalysis_XP_2GQG = pathFolder.analyses("2GQG_XPdock")

# analysis
#dockingAnalysis(psdfDokingXP_2GQG, ltab, pCHEMBLClean, prDockingPoseXP_2HYY, pranalysis_XP_2GQG)

#lprpose = [prDockingPoseXP_2HYY, prDockingPoseXP_3QRJ, prDockingPoseXP_2F4J, prDockingPoseXP_2FO0, prDockingPoseXP_2G2H]




#################
#################
# Analyse NAMS  #
#################
#################

#mcs = MCS.MCSMatrix(ltab, pathFolder.analyses("MCS"))
#mcs.computeMatrixMCS()
#mcs.selectAnalogsMatrix(compoundID="CHEMBL941") # select specificaly a compound

####################
####################
#  Molecular Desc  #
####################
####################

# set variable
CORCOEF = 0.90
prPoses = prDockingPoseXP_2HYY
MAXQUANTILE = 95
Desc1D2D = 1
Desc3D = 1
prDescLigStatic = pathFolder.analyses(psub="DescLig2D3D")

#compute desc and analysis #
############################
#DescLigs = analyzeDBLig.DescriptorLig(ctabAll.table, prPoses, CORCOEF, MAXQUANTILE, Desc1D2D, Desc3D, prDescLigStatic)
#DescLigs.computeDesc()

# generate PNG
#analyzeDBLig.generatePNG()

#analyze descriptor #
#####################
#DescLigs.dendoAffinity("1D2D3D", paff)
#DescLigs.dendoAffinity("1D2D", paff)
#DescLigs.dendoAffinity("3D", paff)



###############
# FPI by pose #
###############
#prFPI = pathFolder.analyses("dockingFPI")
#FPIMatrix(sdocking, pprotein, prFPI)


###############################
# MD based on docking poses   #
###############################

# home
prMD = "/home/aborrel/imitanib/results/MD-ABL/"
#pprotein = "/home/aborrel/imitanib/2hyy_MD.pdb"
#prLig = prDockingPoseSP
pranalysis = "/home/borrela2/imatinib/results/MDanalysis/"

# monster
#prMD = "/data/aborrel/imatinib/results/MD-ABL/"
#pprotein = "/data/aborrel/imatinib/2hyy_MD.pdb"

################
# parameter MD #
################
timeMD = "15000.0"
timeframe = "10.0"
stepWait = 9
nbGPU = 3#maybe integrate in initialization, code clearity
nbCPU = 10
stepFrame = 1# reduce the number of extracted frames
nameLig = "UNK"
water = 0
BSCutoff = 6.0

##########
# Run MD #
##########
pathFolder.createFolder(pranalysis)
#computeMD(prLig, prMD, pprotein, pranalysis, nameLig, BSCutoff, timeMD, timeframe, stepWait, stepFrame, water, nbCPU, nbGPU)


###################
# MD descriptors  #
###################
prMDdesc = "/home/borrela2/imatinib/results/analysis/MDdescriptor/"
pranalysis = "/home/borrela2/imatinib/results/MDanalysis/"
pathFolder.createFolder(prMDdesc)

################################
#   for not bash  to be clean  #
################################
istart = 67
#computeMDdesc(pranalysis, prMDdesc,  istart=istart, iend=istart+1, descLig=1, descBS=0, descFPI=0)


################################
##         for bash           ##
# uncomment for descriptor run #
################################

#import sys
#istart = sys.argv[1]
#if istart[-1] == "&":
#    istart = istart[:-1]
#istart = int(istart)
#computeMDdesc(pranalysis, prMDdesc,  istart=istart, iend=istart+1, descLig=1, descBS=1, descFPI=1)



###################
# Analysis RMSD   #
###################
# RMSD
prMDanalysis = prHome + "imatinib/results/MDanalysis/"
prRMSD = pathFolder.createFolder(prHome + "imatinib/results/analysis/RMSDanalysis/")

#cMD = MD.MD(prMD, pranalysis, water, timeMD, timeframe, stepWait, nbGPU, nbCPU, stepFrame)
#CMDanalysis = MDanalysis.trajectoryAnalysis(cMD, timeMD, timeframe, stepFrame)
#CMDanalysis.histRMSD(paff, prMDanalysis, prRMSD)


# develop the 3D plot  based on density #
#########################################
prdensity = pathFolder.createFolder(prHome + "imatinib/results/analysis/densityMD/")
#CMDanalysis.plotDensity(paff, ["Ki"], 500, prMDdesc, prdensity)


##################################
# Create builder for descriptors #
##################################
prMDdesc = "/home/borrela2/imatinib/results/analysis/MDdescriptor/"
paff ="/home/borrela2/imatinib/results/CHEMBL/AffAllcurated"
pdesc2D = "/home/borrela2/imatinib/results/analysis/DescLig2D3D/1D2D.csv"
pdesc3D = "/home/borrela2/imatinib/results/analysis/DescLig2D3D/3D.csv"
SPLITSET = 0.15



###################################
# develop clustering and quality  #
###################################

CUTOFFACT = 7.5
typeAff = "Ki"
prClustering = pathFolder.analyses("ClusteringDescType")
matrixDescBuilder = builderDescMatrix.Builder(prMDdesc, pdesc2D, pdesc3D, paff, CORCOEF, MAXQUANTILE, SPLITSET, typeAff)

#cclust = clusterDesc.clusterDesc(["Lig2D"], typeAff, CUTOFFACT, prClustering)
#cclust.prep(matrixDescBuilder)
#cclust.clusterize()
#cclust.activityCliff()
#cclust.clusterizeTopActive(50)
#cclust.clusterizeTopActive(20)

#matrixDescBuilder = builderDescMatrix.Builder(prMDdesc, pdesc2D, pdesc3D, paff, CORCOEF, MAXQUANTILE, SPLITSET, typeAff)
#cclust = clusterDesc.clusterDesc(["Lig3D", "Lig2D"], typeAff, CUTOFFACT, prClustering)
#cclust.prep(matrixDescBuilder)
#cclust.clusterize()
#cclust.activityCliff()
#cclust.clusterizeTopActive(50)
#cclust.clusterizeTopActive(20)

#matrixDescBuilder = builderDescMatrix.Builder(prMDdesc, pdesc2D, pdesc3D, paff, CORCOEF, MAXQUANTILE, SPLITSET, typeAff)
#cclust = clusterDesc.clusterDesc(["Lig2D", "Lig3D", "Lig"], typeAff, CUTOFFACT, prClustering)
#cclust.prep(matrixDescBuilder)
#cclust.clusterize()
#cclust.activityCliff()
#cclust.clusterizeTopActive(50)
#cclust.clusterizeTopActive(20)

#matrixDescBuilder = builderDescMatrix.Builder(prMDdesc, pdesc2D, pdesc3D, paff, CORCOEF, MAXQUANTILE, SPLITSET, typeAff)
#cclust = clusterDesc.clusterDesc(["Lig2D", "Lig3D", "Lig", "BS", "FPI"], typeAff, CUTOFFACT, prClustering)
#cclust.prep(matrixDescBuilder)
#cclust.clusterize()
#cclust.activityCliff()
#cclust.clusterizeTopActive(50)
#cclust.clusterizeTopActive(20)


######################
# develop QSAR model #
######################
prMDdesc = "/home/borrela2/imatinib/results/analysis/MDdescriptor/"
paff ="/home/borrela2/imatinib/results/CHEMBL/AffAllcurated"
pdesc2D = "/home/borrela2/imatinib/results/analysis/DescLig2D3D/1D2D.csv"
pdesc3D = "/home/borrela2/imatinib/results/analysis/DescLig2D3D/3D.csv"
#prQSAR = pathFolder.analyses("QSARs")

# settings
varsplit = 0.15


#######################
# only 2D descriptor  #
#######################


####################
### MULTIRUN   #####
####################
#typeAff = "All" #Ki


#QSAR1
#typeAff = "IC50"
#prQSAR = pathFolder.analyses("QSARs1")
#matrixDescBuilder = builderDescMatrix.Builder(prMDdesc, pdesc2D, pdesc3D, paff, CORCOEF, MAXQUANTILE, SPLITSET, typeAff)
#QSARLig2D = QSARModeling.QSARModeling(["Lig2D"], typeAff, prQSAR)
#QSARLig2D.prep(matrixDescBuilder)
#QSARLig2D.runQSARModel()


#typeAff = "Ki" #Ki
#QSARLig2D = QSARModeling.QSARModeling(["Lig2D"], typeAff, prQSAR)
#QSARLig2D.prep(matrixDescBuilder)
#QSARLig2D.runQSARModel()


#QSAR2
#typeAff = "IC50"
#prQSAR = pathFolder.analyses("QSARs2")
#matrixDescBuilder = builderDescMatrix.Builder(prMDdesc, pdesc2D, pdesc3D, paff, CORCOEF, MAXQUANTILE, SPLITSET, typeAff)
#QSARLig2D = QSARModeling.QSARModeling(["Lig2D"], typeAff, prQSAR)
#QSARLig2D.prep(matrixDescBuilder)
#QSARLig2D.runQSARModel()


#typeAff = "Ki" #Ki
#matrixDescBuilder = builderDescMatrix.Builder(prMDdesc, pdesc2D, pdesc3D, paff, CORCOEF, MAXQUANTILE, SPLITSET, typeAff)
#QSARLig2D = QSARModeling.QSARModeling(["Lig2D"], typeAff, prQSAR)
#QSARLig2D.prep(matrixDescBuilder)
#QSARLig2D.runQSARModel()

#QSAR3
#typeAff = "IC50"
#prQSAR = pathFolder.analyses("QSARs3")
#matrixDescBuilder = builderDescMatrix.Builder(prMDdesc, pdesc2D, pdesc3D, paff, CORCOEF, MAXQUANTILE, SPLITSET, typeAff)
#QSARLig2D = QSARModeling.QSARModeling(["Lig2D"], typeAff, prQSAR)
#QSARLig2D.prep(matrixDescBuilder)
#QSARLig2D.runQSARModel()


#typeAff = "Ki" #Ki
#matrixDescBuilder = builderDescMatrix.Builder(prMDdesc, pdesc2D, pdesc3D, paff, CORCOEF, MAXQUANTILE, SPLITSET, typeAff)
#QSARLig2D = QSARModeling.QSARModeling(["Lig2D"], typeAff, prQSAR)
#QSARLig2D.prep(matrixDescBuilder)
#QSARLig2D.runQSARModel()


#QSAR4
#typeAff = "IC50"
#prQSAR = pathFolder.analyses("QSARs4")
#matrixDescBuilder = builderDescMatrix.Builder(prMDdesc, pdesc2D, pdesc3D, paff, CORCOEF, MAXQUANTILE, SPLITSET, typeAff)
#QSARLig2D = QSARModeling.QSARModeling(["Lig2D"], typeAff, prQSAR)
#QSARLig2D.prep(matrixDescBuilder)
#QSARLig2D.runQSARModel()


#typeAff = "Ki" #Ki
#matrixDescBuilder = builderDescMatrix.Builder(prMDdesc, pdesc2D, pdesc3D, paff, CORCOEF, MAXQUANTILE, SPLITSET, typeAff)
#QSARLig2D = QSARModeling.QSARModeling(["Lig2D"], typeAff, prQSAR)
#QSARLig2D.prep(matrixDescBuilder)
#QSARLig2D.runQSARModel()

#QSAR5
#typeAff = "IC50"
#prQSAR = pathFolder.analyses("QSARs5")
#matrixDescBuilder = builderDescMatrix.Builder(prMDdesc, pdesc2D, pdesc3D, paff, CORCOEF, MAXQUANTILE, SPLITSET, typeAff)
#QSARLig2D = QSARModeling.QSARModeling(["Lig2D"], typeAff, prQSAR)
#QSARLig2D.prep(matrixDescBuilder)
#QSARLig2D.runQSARModel()


#typeAff = "Ki" #Ki
#matrixDescBuilder = builderDescMatrix.Builder(prMDdesc, pdesc2D, pdesc3D, paff, CORCOEF, MAXQUANTILE, SPLITSET, typeAff)
#QSARLig2D = QSARModeling.QSARModeling(["Lig2D"], typeAff, prQSAR)
#QSARLig2D.prep(matrixDescBuilder)
#QSARLig2D.runQSARModel()

################
# Regular run  #
################
prQSAR = pathFolder.analyses("QSARs")

###############
### for ki ####
###############
#typeAff = "Ki"

#matrixDescBuilder = builderDescMatrix.Builder(prMDdesc, pdesc2D, pdesc3D, paff, CORCOEF, MAXQUANTILE, SPLITSET, typeAff)
#matrixDescBuilder.dtrain = "/home/borrela2/imatinib/results/analysis/QSARs/Lig2D_Ki/trainSet.csv"
#matrixDescBuilder.dtest = "/home/borrela2/imatinib/results/analysis/QSARs/Lig2D_Ki/testSet.csv"


#2D
#QSARLig2D = QSARModeling.QSARModeling(["Lig2D"], typeAff, prQSAR)
#QSARLig2D.prep(matrixDescBuilder)
#QSARLig2D.runQSARModel()


#2D+3D
#QSARLig2D3D = QSARModeling.QSARModeling(["Lig3D", "Lig2D"], typeAff, prQSAR)
#QSARLig2D3D.prep(matrixDescBuilder)
#QSARLig2D3D.runQSARModel()

#2D+3D+MDlig
#CORCOEF = 0.85
#matrixDescBuilder.corcoef = CORCOEF

#QSARLig2D3D = QSARModeling.QSARModeling(["Lig3D", "Lig2D", "Lig"], typeAff, prQSAR)
#QSARLig2D3D.prep(matrixDescBuilder)
#QSARLig2D3D.runQSARModel()

#2D+3D+MDlig+BS
#CORCOEF = 0.80
#matrixDescBuilder.corcoef = CORCOEF
#QSARLig2D3D = QSARModeling.QSARModeling(["Lig3D", "Lig2D", "Lig", "BS"], typeAff, prQSAR)
#QSARLig2D3D.prep(matrixDescBuilder)
#QSARLig2D3D.runQSARModel()


#2D+3D+MDlig+BS+FPI
#QSARLig2D3D = QSARModeling.QSARModeling(["Lig3D", "Lig2D", "Lig", "BS", "FPI"], typeAff, prQSAR)
#QSARLig2D3D.prep(matrixDescBuilder)
#QSARLig2D3D.runQSARModel()


### for IC50 ####
#################
CORCOEF = 0.90
typeAff = "IC50"
matrixDescBuilder = builderDescMatrix.Builder(prMDdesc, pdesc2D, pdesc3D, paff, CORCOEF, MAXQUANTILE, SPLITSET, typeAff)
matrixDescBuilder.dtrain = "/home/borrela2/imatinib/results/analysis/QSARs/Lig2D_IC50/trainSet.csv"
matrixDescBuilder.dtest = "/home/borrela2/imatinib/results/analysis/QSARs/Lig2D_IC50/testSet.csv"

#2D
#QSARLig2D = QSARModeling.QSARModeling(["Lig2D"], typeAff, prQSAR)
#QSARLig2D.prep(matrixDescBuilder)
#QSARLig2D.runQSARModel()


#2D+3D
QSARLig2D3D = QSARModeling.QSARModeling(["Lig3D", "Lig2D"], typeAff, prQSAR)
QSARLig2D3D.prep(matrixDescBuilder)
QSARLig2D3D.runQSARModel()


#2D+3D+MDlig
QSARLig2D3D = QSARModeling.QSARModeling(["Lig3D", "Lig2D", "Lig"], typeAff, prQSAR)
QSARLig2D3D.prep(matrixDescBuilder)
QSARLig2D3D.runQSARModel()


#2D+3D+MDlig+BS
QSARLig2D3D = QSARModeling.QSARModeling(["Lig3D", "Lig2D", "Lig", "BS"], typeAff, prQSAR)
QSARLig2D3D.prep(matrixDescBuilder)
QSARLig2D3D.runQSARModel()



#2D+3D+MDlig+BS+FPI
QSARLig2D3D = QSARModeling.QSARModeling(["Lig3D", "Lig2D", "Lig", "BS", "FPI"], typeAff, prQSAR)
QSARLig2D3D.prep(matrixDescBuilder)
QSARLig2D3D.runQSARModel()
dd


### for All ####
#################
#matrixDescBuilder = builderDescMatrix.Builder(prMDdesc, pdesc2D, pdesc3D, paff, CORCOEF, MAXQUANTILE, SPLITSET, typeAff)

#typeAff = "All"
#matrixDescBuilder.dtrain = prQSAR + "-".join(["Lig2D"]) + "_" + str(typeAff) + "/trainSetMerged.csv"
#matrixDescBuilder.dtest = prQSAR + "-".join(["Lig2D"]) + "_" + str(typeAff) + "/testSetMerged.csv"

#2D
#QSARLig2D = QSARModeling.QSARModeling(["Lig2D"], typeAff, prQSAR)
#QSARLig2D.prep(matrixDescBuilder)
#QSARLig2D.runQSARModel()

#2D+3D
#QSARLig2D3D = QSARModeling.QSARModeling(["Lig3D", "Lig2D"], typeAff, prQSAR)
#QSARLig2D3D.prep(matrixDescBuilder)
#QSARLig2D3D.runQSARModel()

#2D+3D+MDlig
#QSARLig2D3D = QSARModeling.QSARModeling(["Lig3D", "Lig2D", "Lig"], typeAff, prQSAR)
#QSARLig2D3D.prep(matrixDescBuilder)
#QSARLig2D3D.runQSARModel()

#2D+3D+MDlig+BS
#QSARLig2D3D = QSARModeling.QSARModeling(["Lig3D", "Lig2D", "Lig", "BS"], typeAff, prQSAR)
#QSARLig2D3D.prep(matrixDescBuilder)
#QSARLig2D3D.runQSARModel()

#2D+3D+MDlig+BS+FPI
#QSARLig2D3D = QSARModeling.QSARModeling(["Lig3D", "Lig2D", "Lig", "BS", "FPI"], typeAff, prQSAR)
#QSARLig2D3D.prep(matrixDescBuilder)
#QSARLig2D3D.runQSARModel()



##########################################
# case where we consider the Cell lines  #
##########################################

#pCHEMBL = "/home/aborrel/imitanib/CHEMBL/bioactivity-K562.txt"
#pCHEMBLClean = "/home/aborrel/imitanib/CHEMBL/bioactivity-K562_filtered.txt"


#ltableCpd = CleanCHEMBLFileCellLine(pCHEMBL, pCHEMBLClean)

#mcs = MCS.MCSMatrix(ltableCpd, pathFolder.analyses("MCS-K562"))
#mcs.selectAnalogs(compoundID="CHEMBL941")



###############
# Docking SP  #
###############

#pCHEMBL = "/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862.txt"
#pCHEMBLClean = "/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862_filteredKI.txt"


#ltableCpd = CleanCHEMBLFileProtAff(pCHEMBL, pCHEMBLClean, ["IC50", "Ki", "Kd"], lBAout)

#### docking SP ####
####################


#sdocking = parseSDF.sdf(psdfDoking)
#sdocking.parseSDF()
#sdocking.splitPoses(prDockingPoseSP)
#pdockingAnalysis = pathFolder.analyses("dockingKI")

#dscore = sdocking.get_dockingscore()
#dockingScoreAnalysis(dscore, ltableCpd, pdockingAnalysis)

# specific for a compound
#mcs = MCS.MCSMatrix(ltableCpd, pathFolder.analyses("MCS-K562"))
#mcs.selectAnalogs(compoundID="CHEMBL941")



####################
#### docking XP #### Pb because the sdf file include protein also
####################


#pCHEMBL = "/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862.txt"
#pprotein = "/home/aborrel/imitanib/2hyy_dock.pdb"


# docking parsing #
###################

#sdocking = parseSDF.sdf(psdfDokingXP_2HYY)
#sdocking.parseSDF()
#sdocking.splitPoses(prDockingPoseXP_2HYY)
#dscore = sdocking.get_dockingscore()

# select affinity from CHEMBL

# all affinity #
################
#pCHEMBLClean = "/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862_allAff.txt"
#ltableCpdAll = CleanCHEMBLFileProtAff(pCHEMBL, pCHEMBLClean, ["IC50", "Kd", "Ki"], lBAout)
#pdockingXPAnalysis = pathFolder.analyses("dockingXPAll")
#dockingScoreAnalysis(dscore, ltableCpdAll.table, pCHEMBLClean, pdockingXPAnalysis)

# IC50 #
########
#pCHEMBLClean = "/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862_IC50.txt"
#ltableCpdIC50 = CleanCHEMBLFileProtAff(pCHEMBL, pCHEMBLClean, ["IC50"], lBAout)
#pdockingXPAnalysis = pathFolder.analyses("dockingXPIC50")
#dockingScoreAnalysis(dscore, ltableCpdIC50, pCHEMBLClean, pdockingXPAnalysis)

# Kd #
######
#pCHEMBLClean = "/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862_Kd.txt"
#ltableCpdKd = CleanCHEMBLFileProtAff(pCHEMBL, pCHEMBLClean, ["Kd"], lBAout)
#pdockingXPAnalysis = pathFolder.analyses("dockingXPKd")
#dockingScoreAnalysis(dscore, ltableCpdKd, pCHEMBLClean, pdockingXPAnalysis)

# Ki #
######
#pCHEMBLClean = "/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862_Ki.txt"
#ltableCpdKi = CleanCHEMBLFileProtAff(pCHEMBL, pCHEMBLClean, ["Ki"], lBAout)
#pdockingXPAnalysis = pathFolder.analyses("dockingXPKi")
#dockingScoreAnalysis(dscore, ltableCpdKi, pCHEMBLClean, pdockingXPAnalysis)


#################
#  Clustering   #
#################


# by cluster
#pdesc = pathFolder.analyses("desc") + "CorDesc" + str(CORCOEF) + "/"
#prcluster = pathFolder.analyses("clusterOut")

#mcs = MCS.MCSMatrix(ltab, pathFolder.analyses("MCS"))
#for filein in listdir(pdesc):
#    if search("Table", filein):
        #ccluster = cpdClustering.AnalyseClusterCpd(ltab, pfilecluster=pdesc+filein, proutcluster=prcluster, lprdockingpose=lprpose)
        #mcs.selectCluster(pfilecluster=pdesc+filein, prout=prcluster)#maybe pass in ccluster init -> need to change the folder
        #ccluster.summarize()
        #ccluster.superimposedPoseCluster()
        #ccluster.ShaepMatrix()
        #ccluster.FPIbycluster(pprot=pprotei







