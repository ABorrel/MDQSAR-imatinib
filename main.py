from Tix import _dummyComboBox

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

from os import listdir
from re import search

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



def AnalyseDesc(pdesc, pdata, prout, PCA="1", dendo="1", cormatrix="1", hist="1", corcoef=0.0):

    runExternalSoft.DescAnalysis(pdesc, pdata, prout, corcoef, PCA, cormatrix, hist, dendo)

    return


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

def computeFPIBSBased(cMDs, prout, nameLig):

    dout = {}
    prtempFPI = pathFolder.createFolder(prout + "tempFPI/")
    for nameMD in cMDs.lMD.keys():
        prtempMDFPI = pathFolder.createFolder(prtempFPI + nameMD + "/")

        dout[nameMD] = {}
        prframes = cMDs.lMD[nameMD]["prframe"]

        i = 0
        imax = float(cMDs.MDtime)/float(cMDs.interval)
        #imax = 10 # !!!!!!!!!!!!!!!!!!!!!!!!
        dcFpI = {}
        while i <= imax:
            frameName = "frame_" + str("%05d" % (i)) + ".pdb"
            pframe = prframes + "/" + frameName

            pFPItemp = pathFolder.createFolder(prtempMDFPI + frameName[0:-4] + "/")
            dframe = PDB.PDB(pframe, hydrogen=1)

            CFPI = FPI.ligFPI(dframe, pFPItemp, ligID=nameLig)
            CFPI.computeFPI()
            dcFpI[frameName] = CFPI

            i += 1

        FPIMD = FPI.CompareFPIMD(dcFpI, prtempMDFPI)
        FPIMD.MDprop()
        FPIMD.pobaFPI()
        #FPIMD.MDtanimoto() # useless if only ligand is considered


    return dout


##########
#  MAIN  #
##########


###############
#  CONSTANT   #
###############

lBAout = ["CHEMBL3705971"]
lBAout = []

# gleevec = CHEMBL941
# outlier = CHEMBL2382016


###################
# Docking folder  #
###################

pprotein = "/home/aborrel/imitanib/2hyy_dock.pdb"

# SP #
######
psdfDokingSP = "/home/aborrel/imitanib/docking/dockingSP_2hyy/dockingpose.sdf"
prDockingPoseSP = "/home/aborrel/imitanib/results/dockingposeSP/"

# XP #
######
psdfDokingXP = "/home/aborrel/imitanib/docking/dockingXP_2hyy/PoseXP.sdf"
prDockingPoseXP = "/home/aborrel/imitanib/results/dockingposeXP/"


# case where we consider the binding affinity #
###############################################

pCHEMBL = "/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862.txt"
pCHEMBLClean = "/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862_filtered.txt"


# Load table from CHEMBL #
##########################
ltab = CleanCHEMBLFileProtAff(pCHEMBL, pCHEMBLClean, ["IC50", "Ki", "Kd"], lBAout)


#################
# Analyse NAMS  #
#################
#mcs = MCS.MCSMatrix(ltab, pathFolder.analyses("MCS"))
#mcs.computeMatrixMCS()
#mcs.selectAnalogsMatrix(compoundID="CHEMBL941") # select specificaly a compound


####################
#  Molecular Desc  #
####################
pdesc = pathFolder.analyses(psub="desc") + "tableDesc"
plog = pathFolder.analyses(psub="desc") + "log.txt"
pdescglobal = liganddescriptors.MolecularDesc(ltab, pdesc, prDockingPoseXP, plog)# add docking XP
#AnalyseDesc(pdesc, pCHEMBLClean, pathFolder.analyses("desc"), corcoef=0.7)




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

# monster
#prMD = "/data/aborrel/imatinib/results/MD-ABL/"
#pprotein = "/data/aborrel/imatinib/2hyy_MD.pdb"


# parameter MD
timeMD = "15000.0"
timeframe = "10.0"
stepWait = 16

# 1. Merge poses and proteins
#cMDs = MD.MD(prMD, timeMD, timeframe, stepWait)
#cMDs.initialisation(prDockingPoseSP, pprotein)
#cMDs.runMultipleMD() # run MD

# 2. analyse MD
# name ligand for the MD
namelig = "UNK"# classic name given by glide

#cMDs.extractFrame()
#cMDs.analyseAllMD(RMSD=1, ligAnalysis=0, nameLig=namelig)


# 3. FPI computation
# ligand + BS based
#prFPI = pathFolder.analyses("MD_FPI")
#computeFPIBSBased(cMDs, prFPI, namelig)





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
#pdesc = pathFolder.analyses("desc")
#prcluster = pathFolder.analyses("clusterOut")

#mcs = MCS.MCSMatrix(ltableCpd, pathFolder.analyses("MCS"))
#for filein in listdir(pdesc):
#    if search("Table", filein):
#        ccluster = cpdClustering.AnalyseClusterCpd(pfilecluster=pdesc+filein, proutcluster=prcluster, prdockingpose=prDockingPoseXP)
#        ccluster.superimposedPoseCluster()
#        mcs.selectCluster(pfilecluster=pdesc+filein, prout=prcluster)#maybe pass in ccluster init -> need to change the folder
#        ccluster.ShaepMatrix()
#        ccluster.FPIbycluster(pprot=pprotein)

