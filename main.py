import tableParse
import liganddescriptors
import pathFolder
import runExternalSoft



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

##########
#  MAIN  #
##########

pCHEMBL = "/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862.txt"
pCHEMBLClean = "/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862_filtered.txt"

#ltab = CleanCHEMBLFile(pCHEMBL, pCHEMBLClean)

# matrix of MCS #




pdesc = pathFolder.analyses(psub="desc") + "tableDesc.csv"
plog = pathFolder.analyses(psub="desc") + "log.txt"
#MolecularDesc(ltab, pdesc, plog)

AnalyseDesc(pdesc, pCHEMBLClean, pathFolder.analyses("desc"), corcoef=0.7)

psdfDoking = "/home/aborrel/imitanib/results/dockingpose.sdf"

