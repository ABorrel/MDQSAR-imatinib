import tableParse
import liganddescriptors
import pathFolder
import runExternalSoft
import MCS
import parseSDF
import FPI
import PDB


def CleanCHEMBLFile(pfilin, pfilout):

    # add short cut if filtered table exist !!!!!

    table = tableParse.CHEMBL(pfilin)
    table.parseCHEMBLFile()
    print len(table.table), "Init"

    table.selectConfidencecore(cutoff=9)
    print len(table.table), "prot confidence"

    table.getOnlyExactConstant()
    print len(table.table), "strict value"

    #table.getOnlyIC50()
    #print len(table.table), "IC50"

    table.MergeIdenticCHEMBLIDforACtivity()
    print len(table.table), "Repetition"

    table.writeTable(pfilout)


    return table.table


def MolecularDesc(ltable, pfilout, plog):

    logfile = open(plog, "w")

    for compound in ltable:
        dcompound = liganddescriptors.Descriptors(compound, logfile)

        if dcompound.log == "ERROR":
            continue

        dcompound.get_descriptorOD1D()
        dcompound.get_descriptor2D()

        dcompound.writeTablesDesc(pfilout)


def AnalyseDesc(pdesc, pdata, prout, PCA="1", dendo="1", cormatrix="1", hist="1", corcoef=0.0):

    runExternalSoft.DescAnalysis(pdesc, pdata, prout, corcoef, PCA, cormatrix, hist, dendo )

    return


def dockingScoreAnalysis(ddockingscore, ltabCHEMBL, prout):

    pfilout = prout + "ScoreVSAff.txt"
    filout = open(pfilout, "w")
    filout.write("IDCHEMBL\tDock_score\tAff\n")

    for daff in ltabCHEMBL:
        try: filout.write(str(daff["CMPD_CHEMBLID"]) + "\t" + str(ddockingscore[daff["CMPD_CHEMBLID"]]) + "\t" + str(daff["PCHEMBL_VALUE"]) + "\n")
        except: pass
    filout.close()

    runExternalSoft.corPlot(pfilout)


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



##########
#  MAIN  #
##########

pCHEMBL = "/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862.txt"
pCHEMBLClean = "/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862_filtered.txt"

ltab = CleanCHEMBLFile(pCHEMBL, pCHEMBLClean)

# matrix of MCS #
#mcs = MCS.MCSMatrix(ltab)
#mcs.computeMatrixMCS(pathFolder.analyses("MCS"))



#pdesc = pathFolder.analyses(psub="desc") + "tableDesc.csv"
#plog = pathFolder.analyses(psub="desc") + "log.txt"
#MolecularDesc(ltab, pdesc, plog)

#AnalyseDesc(pdesc, pCHEMBLClean, pathFolder.analyses("desc"), corcoef=0.7)

psdfDoking = "/home/aborrel/imitanib/results/dockingpose.sdf"
prDockingPose = "/home/aborrel/imitanib/results/dockingpose/"
prFPI = pathFolder.analyses("dockingFPI")
pdockingAnalysis = pathFolder.analyses("docking")
pprotein = "/home/aborrel/imitanib/2hyy_dock.pdb"

sdocking = parseSDF.sdf(psdfDoking)
sdocking.parseSDF()
sdocking.splitPoses(prDockingPose)
dscore = sdocking.get_dockingscore()

dockingScoreAnalysis(dscore, ltab, pdockingAnalysis)

# FPI by pose
FPIMatrix(sdocking, pprotein, prFPI)

