import tableParse
import liganddescriptors
import pathFolder
import runExternalSoft
import MCS
import parseSDF
import MD
import FPI
import PDB
import cpdClustering
import MDdescriptors
import MDQSAR

from os import listdir, makedirs
from re import search


#########################
##   MAIN FUNCTIONS    ##
#########################


def CleanCHEMBLFileProtAff(pfilin, pfilout, ltypeAff, lBAout):

    # add short cut if filtered table exist !!!!!

    table = tableParse.CHEMBL(pfilin)
    table.parseCHEMBLFile()
    print len(table.table), "Init"

    table.selectConfidencecore(cutoff=9)
    print len(table.table), "prot confidence"

    table.getOnlyExactConstant()
    print len(table.table), "strict value"

    table.getByTypeOfAff(ltypeAff)
    print len(table.table), ltypeAff

    table.MergeIdenticCHEMBLIDforACtivity()
    print len(table.table), "Repetition"

    table.selectAssayType("B")
    print len(table.table), "Type assay"

    # remove some biassay
    table.removeBA(lBAout)
    print len(table.table)

    table.writeTable(pfilout)


    return table.table



def CleanCHEMBLFileCellLine(pfilin, pfilout, ltypeaff=["IC50"]):

    # add short cut if filtered table exist !!!!!

    table = tableParse.CHEMBL(pfilin)
    table.parseCHEMBLFile()
    print len(table.table), "Init"

    table.getOnlyExactConstant()
    print len(table.table), "strict value"

    table.selectAssayType("F")
    print len(table.table), "Type assay"

    table.MergeIdenticCHEMBLIDforACtivity()
    print len(table.table), "Repetition"

    table.getByTypeOfAff(ltypeaff)
    print len(table.table), ltypeaff

    table.writeTable(pfilout)


    return table.table



def AnalyseDesc(pdesc, pdata, prout, PCA="1", dendo="1", cormatrix="1", hist="1", clustering="1", corcoef=0.0):

    prout = prout + "CorDesc" + str(corcoef) + "/"
    try: makedirs(prout)
    except:pass
    runExternalSoft.DescAnalysis(pdesc, pdata, prout, corcoef, PCA, cormatrix, hist, dendo, clustering)

    return prout


def dockingScoreAnalysis(ddockingscore, ltabCHEMBL, ptableCHEMBL, prout):

    pfilout = prout + "ScoreVSAff.txt"
    filout = open(pfilout, "w")
    filout.write("IDCHEMBL\tDock_score\temodel\tAff\n")

    for daff in ltabCHEMBL:
        try: filout.write(str(daff["CMPD_CHEMBLID"]) + "\t" + str(ddockingscore[daff["CMPD_CHEMBLID"]]["r_i_docking_score"])
                          + "\t" + str(ddockingscore[daff["CMPD_CHEMBLID"]]["r_i_glide_emodel"])
                          + "\t" + str(daff["PCHEMBL_VALUE"]) + "\n")
        except: pass
    filout.close()

    runExternalSoft.corPlot(pfilout, ptableCHEMBL)


def FPIMatrix(sdocking, pprotein, prFPI):

    pmatrixFPI = prFPI + "MFPI.txt"

    i = 0
    imax = len(sdocking.lposefiles)
    imax = 3
    while i < imax:
        j = i + 1
        while j < imax:
            cprot = PDB.PDB(PDB_input=pprotein, hydrogen=1)

            pligPDBi = runExternalSoft.babelConvertSDFtoPDB(sdocking.lposefiles[i])
            pligPDBj = runExternalSoft.babelConvertSDFtoPDB(sdocking.lposefiles[j])

            print pligPDBi
            cposei = PDB.PDB(PDB_input=pligPDBi, hydrogen=1)
            cposej = PDB.PDB(PDB_input=pligPDBj, hydrogen=1)

            sFPIi = FPI.ligFPI(cPDB=cprot, ligin=cposei, prFPI=prFPI)
            sFPIj = FPI.ligFPI(cPDB=cprot, ligin=cposej, prFPI=prFPI)

            dout = sFPIi.compareFPI(sFPIj)

            print dout

            j = j + 1
        i = i + 1


def dockingAnalysis(psdfDoking, ltableCpd, ptableCpd, prpose, pranalysis):

    sdocking = parseSDF.sdf(psdfDoking)
    sdocking.parseSDF()
    sdocking.splitPoses(prpose)

    dscore = sdocking.get_dockingscore()
    dockingScoreAnalysis(dscore, ltableCpd, ptableCpd, pranalysis)
    return dscore



def computeMD(prLigand, prMD, pprotein, nameLig, BSCutoff, timeMD, timeframe, stepWait, stepFrame, water, nbCPU, nbGPU):
    """

    :param prLigand: folder of ligands or poses
    :param pprotein: protein in PDB already prepared using the maestro
    :param timeMD: time MD in ps
    :param timeframe: time between two frame
    :param stepWait: time to wait for multiprocess
    :param stepFrame: frame extraction
    :param water: extract water molecules
    :param nbCPU: CPU max in parrallel
    :param nbGPU: GPU max in paralele, random selection
    :return:
    """
    # 1. Merge poses and proteins
    cMDs = MD.MD(prMD, pranalysis, water, timeMD, timeframe, stepWait, nbGPU, nbCPU, stepFrame)
    cMDs.initialisation(prLigand, pprotein)
    cMDs.runMultipleMD()# run MD

    # 2. Extract frames
    # extract frame
    cMDs.centerFrame()
    cMDs.extractFrame()

    # extract BS and ligand
    cMDs.extractLigBSbyFrame(BSCutoff, nameLig, clean=0)

    # 3. compute RMSD
    cMDs.analyseRMSD()






def computeMDdesc(prMD, prout):

    lpMDresult = listdir(prMD)

    dcMDdesc = {}
    for MDresult in lpMDresult:
        jobname = MDresult.split("_")[0]
        prlig = prMD + MDresult + "/lig/"
        prBS = prMD + MDresult  + "/BSs/"
        prframes = prMD + MDresult + "/framesMD/"
        prMDdesc = prout

        # compute different descriptors
        cMD = MDdescriptors.MDdescriptors(jobname, prlig, prBS, prframes, prMDdesc)
        cMD.computeLigDesc()
        cMD.computeBSDesc()
        cMD.computeFPI()

        dcMDdesc[jobname] = cMD




def computeMDQSAR(prDesc, ltypedesc, prout):

    cMDQSAR = MDQSAR.MDQSAR(prDesc, prout)
    cMDQSAR.builtDataset()

    cMDQSAR.runQSARModel(ltypedesc)

    return




##########
#  MAIN  #
##########

# case where we consider the binding affinity #
###############################################

###############
#  CONSTANT   #
###############

lBAout = ["CHEMBL3705971"]
lBAout = []

CORCOEF = 0.8
# gleevec = CHEMBL941
# outlier = CHEMBL2382016


################
# TABLE CHEMBL #
################

pCHEMBL = "/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862.txt"
pCHEMBLClean = "/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862_filtered.txt"
#ltab = CleanCHEMBLFileProtAff(pCHEMBL, pCHEMBLClean, ["IC50", "Ki", "Kd"], lBAout)



###################
# Docking folder  #
###################

pprotein = "/home/aborrel/imitanib/2hyy_dock.pdb"

# SP-2HYY #
###########
psdfDokingSP = "/home/aborrel/imitanib/docking/dockingSP_2hyy/dockingpose.sdf"
prDockingPoseSP = "/home/aborrel/imitanib/results/dockingposeSP_2HYY/"
pprotein_2HYY = "/home/aborrel/imitanib/protein/2HYY_dock.pdb"
pranalysis_SP_2HYY = pathFolder.analyses("2HYY_SPdock")

# analysis
#dockingAnalysis(psdfDokingSP, ltab, pCHEMBLClean, prDockingPoseSP, pranalysis_SP_2HYY)


# XP #
######
# native - 2HYY #
psdfDokingXP_2HYY = "/home/aborrel/imitanib/docking/dockingXP_2hyy/PoseXP.sdf"
prDockingPoseXP_2HYY = "/home/aborrel/imitanib/results/dockingposeXP_2HYY/"
pprotein_2HYY = "/home/aborrel/imitanib/protein/2HYY_dock.pdb"
pranalysis_XP_2HYY = pathFolder.analyses("2HYY_XPdock")

# analysis
#dockingAnalysis(psdfDokingXP_2HYY, ltab, pCHEMBLClean, prDockingPoseXP_2HYY, pranalysis_XP_2HYY)


# mutated 3QRJ
psdfDokingXP_3QRJ = "/home/aborrel/imitanib/docking/dockingXP_3QRJ/3QRJ_XP_pose.sdf"
prDockingPoseXP_3QRJ = "/home/aborrel/imitanib/results/dockingposeXP_3QRJ/"
pprotein_3QRJ = "/home/aborrel/imitanib/protein/3QRJ_dock.pdb"
pranalysis_XP_3QRJ = pathFolder.analyses("3QRJ_XPdock")

# analysis
#dockingAnalysis(psdfDokingXP_3QRJ, ltab, pCHEMBLClean, prDockingPoseXP_3QRJ, pranalysis_XP_3QRJ)


# mutated 2F4J
psdfDokingXP_2F4J = "/home/aborrel/imitanib/docking/dockingXP_2F4J/2F4J_XP_pose.sdf"
prDockingPoseXP_2F4J = "/home/aborrel/imitanib/results/dockingposeXP_2F4J/"
pprotein_2F4J = "/home/aborrel/imitanib/protein/2F4J_dock.pdb"
pranalysis_XP_2F4J = pathFolder.analyses("2F4J_XPdock")

# analysis
#dockingAnalysis(psdfDokingXP_2F4J, ltab, pCHEMBLClean, prDockingPoseXP_2F4J, pranalysis_XP_2F4J)


# mutated 2FO0
psdfDokingXP_2FO0 = "/home/aborrel/imitanib/docking/dockingXP_2FO0/2FO0_XP_pose.sdf"
prDockingPoseXP_2FO0 = "/home/aborrel/imitanib/results/dockingposeXP_2FO0/"
pprotein_2FO0 = "/home/aborrel/imitanib/protein/2FO0_dock.pdb"
pranalysis_XP_2FO0 = pathFolder.analyses("2FO0_XPdock")

# analysis
#dockingAnalysis(psdfDokingXP_2FO0, ltab, pCHEMBLClean, prDockingPoseXP_2FO0, pranalysis_XP_2FO0)

# mutated 2G2H
psdfDokingXP_2G2H = "/home/aborrel/imitanib/docking/dockingXP_2G2H/2G2H_XP_pose.sdf"
prDockingPoseXP_2G2H = "/home/aborrel/imitanib/results/dockingposeXP_2G2H/"
pprotein_2G2H = "/home/aborrel/imitanib/protein/2G2H_dock.pdb"
pranalysis_XP_2G2H = pathFolder.analyses("2G2H_XPdock")

# analysis
#dockingAnalysis(psdfDokingXP_2G2H, ltab, pCHEMBLClean, prDockingPoseXP_2G2H, pranalysis_XP_2G2H)


# mutated 2GQG - verif
#psdfDokingXP_2GQG = "/home/aborrel/imitanib/docking/dockingXP_2GQG/2GQG_XP_pose.sdf"
#prDockingPoseXP_2GQG = "/home/aborrel/imitanib/results/dockingposeXP_2GQG/"
#pprotein_2GQG = "/home/aborrel/imitanib/protein/2GQG_dock.pdb"
#pranalysis_XP_2GQG = pathFolder.analyses("2GQG_XPdock")

# analysis
#dockingAnalysis(psdfDokingXP_2GQG, ltab, pCHEMBLClean, prDockingPoseXP_2HYY, pranalysis_XP_2GQG)

lprpose = [prDockingPoseXP_2HYY, prDockingPoseXP_3QRJ, prDockingPoseXP_2F4J, prDockingPoseXP_2FO0, prDockingPoseXP_2G2H]


#################
# Analyse NAMS  #
#################
#mcs = MCS.MCSMatrix(ltab, pathFolder.analyses("MCS"))
#mcs.computeMatrixMCS()
#mcs.selectAnalogsMatrix(compoundID="CHEMBL941") # select specificaly a compound


####################
#  Molecular Desc  #
####################
#pdesc = pathFolder.analyses(psub="desc") + "tableDesc"
#plog = pathFolder.analyses(psub="desc") + "log.txt"
#pdescglobal = liganddescriptors.MolecularDesc(ltab, pdesc, prDockingPoseXP, plog)# add docking XP
#AnalyseDesc(pdescglobal, pCHEMBLClean, pathFolder.analyses("desc"), corcoef=CORCOEF)




##################################
#  Docking analysis -SP docking  #
##################################


# monster #
###########
#pprotein = "/data/aborrel/imatinib/2hyy_dock.pdb"
#psdfDoking = "/data/aborrel/imatinib/results/dockingpose.sdf"
#prDockingPoseSP = "/data/aborrel/imatinib/results/dockingposeSP/"

#sdocking = parseSDF.sdf(psdfDoking)
#sdocking.parseSDF()
#sdocking.splitPoses(prDockingPoseSP)
#pdockingAnalysis = pathFolder.analyses("docking")


#dscore = sdocking.get_dockingscore()
#dockingScoreAnalysis(dscore, ltab, pdockingAnalysis)



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
pprotein = "/home/aborrel/imitanib/2hyy_MD.pdb"
prLig = prDockingPoseSP


# monster
#prMD = "/data/aborrel/imatinib/results/MD-ABL/"
#pprotein = "/data/aborrel/imatinib/2hyy_MD.pdb"


# parameter MD #
################
timeMD = "15000.0"
timeframe = "10.0"
stepWait = 9
nbGPU = 3#maybe integrate in initialization, code clearity
nbCPU = 10
stepFrame = 10# reduce the number of extracted frames
nameLig = "UNK"
water = 0
BSCutoff = 6.0

# Run MD #
##########
#computeMD(prLig, prMD, pprotein, nameLig, BSCutoff, timeMD, timeframe, stepWait, stepFrame, water, nbCPU, nbGPU)




###################
# MD descriptors  #
###################
prMDdesc = "/home/aborrel/imitanib/results/analysis/MD_descriptor/"
pranalysis = "/home/aborrel/imitanib/results/MDanalysis/"
pathFolder.createFolder(prMDdesc)
#computeMDdesc(pranalysis, prMDdesc)





######################
# develop QSAR model #
######################

prQSAR = pathFolder.analyses("MDQSARs")
computeMDQSAR(prMDdesc, ["Lig"], prQSAR)







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


pCHEMBL = "/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862.txt"
pprotein = "/home/aborrel/imitanib/2hyy_dock.pdb"


# docking parsing #
###################

#sdocking = parseSDF.sdf(psdfDoking)
#sdocking.parseSDF()
#sdocking.splitPoses(prDockingPoseXP)
#dscore = sdocking.get_dockingscore()

# select affinity from CHEMBL

# all affinity #
################
#pCHEMBLClean = "/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862_allAff.txt"
#ltableCpdAll = CleanCHEMBLFileProtAff(pCHEMBL, pCHEMBLClean, ["IC50", "Kd", "Ki"], lBAout)
#pdockingXPAnalysis = pathFolder.analyses("dockingXPAll")
#dockingScoreAnalysis(dscore, ltableCpdAll, pCHEMBLClean, pdockingXPAnalysis)

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
        #ccluster.FPIbycluster(pprot=pprotein)

