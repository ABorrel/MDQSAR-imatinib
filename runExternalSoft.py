from os import system, path, remove, getcwd, chdir
from shutil import copyfile
from re import search
import toolbox
from time import sleep

LIGPREP = "/opt/schrodinger2016-4/ligprep"
PADEL = "/home/aborrel/softwares/padel/PaDEL-Descriptor.jar"
STRUCTCONVERT = "/opt/schrodinger2016-4/utilities/structconvert"
MULTISMIM = "/opt/schrodinger2016-4/utilities/multisim"
STRUCTCAT = "/opt/schrodinger2016-4/utilities/structcat"


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


def convertPDBtoMAE(pPDB):

    pMAE = pPDB[:-4] + ".mae"
    cmdconvert = STRUCTCONVERT + " -ipdb " + str(pPDB) + " -omae " + pMAE
    print cmdconvert
    system(cmdconvert)

    if path.exists(pMAE):
        return pMAE
    else:
        print "ERROR---"
        return "ERROR"


def multisimSystemBuilder(jobname, pMAE, WAIT=1):

    pcms = pMAE[:-4] + ".cms"
    # cp script
    pmsj = path.dirname(pMAE) + "/systembuilder.msj"
    copyfile("./systembuilder.msj", pmsj)

    # move in directory
    prsource = getcwd()
    prrun = path.dirname(pMAE)
    chdir(prrun)
    print getcwd()

    if WAIT == 1 :
        cmdMUL = MULTISMIM + " -JOBNAME " + str(jobname) + " -m " + str(pmsj) + " " + str(pMAE) + " -o " + pcms + " -WAIT"
        print cmdMUL
        system(cmdMUL)
    else:
        cmdMUL = MULTISMIM + " -JOBNAME " + str(jobname) + " -m " + str(pmsj) + " " + str(pMAE) + " -o " + pcms
        " -WAIT"
        print cmdMUL
        system(cmdMUL)

    if not path.exists(pcms):
        # go back in source folder
        chdir(prsource)
        return "ERROR - System Builder"
    else:
        # go back in source folder
        chdir(prsource)
        return pcms


def multisimGDesmond(jobname, pcms, timeMDns, frameInverval, WAIT=1):

    prDM = path.dirname(pcms) + "/"
    pmsj = prDM + jobname + ".msj"
    pcfg = prDM + jobname + ".cfg"

    # write file
    toolbox.writeFilesParamaterDesmond(pmsj, pcfg, timeMDns, frameInverval)

    # move on run folder
    prsource = getcwd()
    prrun = path.dirname(pcms)
    chdir(prrun)


    cmdDesmond = MULTISMIM + " -JOBNAME " + str(jobname) + " -maxjob 1 -cpu 1 -m " + str(pmsj) + " -c " + str(pcfg) + " " + str(pcms) + " -mode umbrella -set stage[1].set_family.md.jlaunch_opt=[\"-gpu\"] -o " + pcms[:-4] + "-out.cms"

    if WAIT == 1:
        cmdDesmond = cmdDesmond + " -WAIT"
    print cmdDesmond
    system(cmdDesmond)
    chdir(prsource)


def concateneStructure(pPDBprot, pligSDF, pmaeComplex):

    cmdStructcat = STRUCTCAT + " " + str(pPDBprot) + " " + str(pligSDF) + " -omae " + str(pmaeComplex)
    print cmdStructcat
    system(cmdStructcat)

    if path.exists(pmaeComplex) and path.getsize(pmaeComplex) > 0:
        return pmaeComplex
    else:
        dddd
        return "ERROR"





