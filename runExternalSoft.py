from os import system, path, remove
from re import search
from time import sleep

LIGPREP = "/opt/schrodinger2016-4/ligprep"
PADEL = "/home/aborrel/softwares/padel/PaDEL-Descriptor.jar"


def runLigprep(psmilin, forcefield="OPLS3", stereoisoster=1):

    """Maybe fix more option"""

    if forcefield == "OPLS3":
        bff = "16"
    else:
        bff = "14"

    cmd = LIGPREP + " -ismi " + psmilin + " -osd " + psmilin[0:-4] + ".sdf" + " -bff " + str(bff) + " -epik -s " + str(stereoisoster) + " -WAIT -NJOBS 3"

    print cmd
    system(cmd)

    # control if file exist
    if not path.exists(psmilin[0:-4] + ".sdf"):
        return "Ligprep ERROR"
    else:
        try:remove("tem.log")
        except:pass
        try:remove("tem-dropped.smi")
        except: pass
        try:remove("tem-dropped-indices.txt")
        except: pass

    return psmilin[0:-4] + ".sdf"

def runPadel(prin=""):
    """Input include a folder of sdf file"""
    if prin == "":
        return "ERROR - Padel Input"
    else:
        cmd = "java -jar " + PADEL + " -maxruntime 100 -3d -dir " + str(prin) + " -file " + prin + "tem.desc"
        print cmd
        system(cmd)

    return prin + "tem.desc"


def babelConvertSDFtoSMILE(sdfread, clean_smi=0, rm_smi=1):

    tempsdf = open("tempsdf.sdf", "w")
    tempsdf.write(sdfread)
    tempsdf.close()

    psmile = "tempsmile.smi"

    cmd_convert = "babel tempsdf.sdf " + psmile + " 2>/dev/null"
    system(cmd_convert)

    try : filin = open (psmile, "r")
    except : return "0"
    l_Fline = filin.readlines()
    filin.close()
    try : smile = l_Fline[0].split("\t")[0]
    except : return "0"

    # rewrite path in filout
    if clean_smi == 1:
        filout = open (psmile, "w")
        filout.write (str (smile))
        filout.close ()

    if rm_smi == 1:
        system("rm " + psmile)


    return smile



def PCAplot(pfildesc, pfildata, corcoef, prout):

    cmdplotPCA = "./PCAplot.R " + str(pfildesc) + " " + str(pfildata) + " " + str(corcoef) + " " + str(prout)

    print cmdplotPCA
    system(cmdplotPCA)

    return

def DescAnalysis(pdesc, paffinity, prout, valcor, PCA, corMatrix, hist, dendo):

    cmdVisu = "./visualization.R " + str(pdesc) + " " + str(paffinity) + " " + str(prout) + " " + str(valcor) + " " + \
        str(PCA) + " " + str(corMatrix) + " " + str(hist) + " " + str(dendo)

    print cmdVisu
    system(cmdVisu)


def MatrixMCS(pmatrix, paff, ptext):

    cmdmatrix = "./matrixMCS.R " + pmatrix + " " + paff + " " + ptext
    print cmdmatrix
    system(cmdmatrix)


def MDSMCS(pmatrix, paff):

    cmdMDS = "./MDSMCS.R " + pmatrix + " " + paff
    print cmdMDS
    system(cmdMDS)


def corPlot(pfilin):

    cmdCor = "./corplot.R " + pfilin
    print cmdCor
    system(cmdCor)


def babelConvertSDFtoPDB(psdf):

    cmdconvert = "babel " + psdf + " " + psdf[:-4] + ".pdb 2>/dev/null"
    print cmdconvert
    system(cmdconvert)

    return psdf[:-4] + ".pdb"
