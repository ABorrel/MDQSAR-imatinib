from os import system, path, remove, getcwd, chdir, makedirs, listdir, removedirs
from shutil import copyfile
from re import search
from time import sleep
import subprocess
import sys


import toolbox
import pathFolder
# import parser
sys.path.insert(0, "./Parser/") # for window dev
import PDB


LIGPREP = "/opt/schrodinger2017-1/ligprep"
PADEL = "/home/aborrel/softwares/padel/PaDEL-Descriptor.jar"
TMalign = "/home/aborrel/softwares/TMalign/TMalign"
SHAEP = "/home/aborrel/softwares/shaep/shaep"
PRNACCESS = "/home/aborrel/softwares/NACCESS/"
PRADI = "/home/aborrel/softwares/RADI.4.0.1"

# for home computer
STRUCTCONVERT = "/opt/schrodinger2017-3/utilities/structconvert"
MULTISMIM = "/opt/schrodinger2017-3/utilities/multisim"
STRUCTCAT = "/opt/schrodinger2017-3/utilities/structcat"
RUN = "/opt/schrodinger2017-3/run"
KRAKENX = "/home/aborrel/softwares/KRAKENX/dist/KrakenX.jar"

# for monster
#LIGPREP = "/opt/schrodinger2016-3/ligprep"
#STRUCTCONVERT = "/opt/schrodinger2016-3/utilities/structconvert"
#MULTISMIM = "/opt/schrodinger2016-3/utilities/multisim"
#STRUCTCAT = "/opt/schrodinger2016-3/utilities/structcat"
#RUN = "/opt/schrodinger2016-3/run"


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


def runPadel(prin, pfilout, force=0):
    """Input include a folder of sdf file"""
    if force == 1:
        try:
            remove(pfilout)
        except:pass

    if path.exists(pfilout):
        filout = open(pfilout, "r")
        llines = filout.readlines()
        filout.close()
        if len(listdir(prin)) <= len(llines):
            return pfilout

    if not path.exists(prin):
        return "ERROR - Padel Input"
    else:
        cmd = "java -jar " + PADEL + " -maxruntime 5000 -3d -dir " + str(prin) + " -file " + pfilout
        print cmd
        system(cmd)

    return pfilout


# too slow !!!
def runKrakenX(prtemp, lcpd, pdesc):

    if path.exists(pdesc) and path.getsize(pdesc) > 100:
        return pdesc

    else:

        filparmIn = open("params_template.txt", "r")
        parmKrakenX = filparmIn.read()
        filparmIn.close()

        lst = open(prtemp + "sdf.lst", "w")
        lst.write("\n".join(lcpd))
        lst.close()

        pparm = prtemp + "param.txt"
        filparTmp = open(pparm, "w")
        filparTmp.write(parmKrakenX%(prtemp + "sdf.lst", pdesc))
        filparTmp.close()
        cmdKrakenX = "java -jar " + KRAKENX + " " + pparm

        prsource = getcwd()
        chdir(prtemp)
        print "RUN ---- "
        print cmdKrakenX
        system(cmdKrakenX)
        chdir(prsource)

    #return pdesc


def babelConverttoSDF(ppdb, pfilout =""):

    if pfilout == "":
        pfilout = ppdb[:-4] + ".sdf"

    cmdconvert = "/usr/bin/babel " + ppdb + " " +  pfilout + " 2>/dev/null"
    system(cmdconvert)

    return pfilout



def babelConvertSDFtoSMILE(sdfread, clean_smi=0, rm_smi=1):

    tempsdf = open("tempsdf.sdf", "w")
    tempsdf.write(sdfread)
    tempsdf.close()

    psmile = "tempsmile.smi"

    cmd_convert = "babel tempsdf.sdf " + psmile + " 2>/dev/null"
    system(cmd_convert)

    try: filin = open (psmile, "r")
    except: return "0"
    l_Fline = filin.readlines()
    filin.close()
    try: smile = l_Fline[0].split("\t")[0]
    except: return "0"

    # rewrite path in filout
    if clean_smi == 1:
        filout = open(psmile, "w")
        filout.write(str(smile))
        filout.close()

    if rm_smi == 1:
        system("rm " + psmile)


    return smile



def PCAplot(pfildesc, pfildata, corcoef, prout):

    cmd = "./PCAplot.R " + str(pfildesc) + " " + str(pfildata) + " " + str(corcoef) + " " + str(prout)
    runRscript(cmd)

    return

def DescAnalysis(pdesc, paffinity, prout, valcor, maxQuantile, PCA, corMatrix, hist, dendo, clustering):

    cmdVisu = "./visualization.R " + str(pdesc) + " " + str(paffinity) + " " + str(prout) + " " + str(valcor) + " " + \
         str(maxQuantile) + " " + str(PCA) + " " + str(corMatrix) + " " + str(hist) + " " + str(dendo) + " " + str(clustering)

    runRscript(cmdVisu)


def MatrixMCS(pmatrix, paff, ptext, probmatrix=1):

    cmdmatrix = "./matrixMCS.R " + pmatrix + " " + paff + " " + ptext + " " + str(probmatrix)
    runRscript(cmdmatrix)


def MDSMCS(pmatrix, paff):

    cmdMDS = "./MDSMCS.R " + pmatrix + " " + paff
    runRscript(cmdMDS)


def corPlot(pfilin, pchembl, prout, typeplot="dockAff"):

    cmdCor = "./corplot.R " + pfilin + " " + pchembl + " " + prout + " " + typeplot
    runRscript(cmdCor)


def histAffinity(paff):

    cmdHist = "./distributionAff.R " + str(paff)
    runRscript(cmdHist)


def histRMSD(pRMSD, prout):

    cmdHist = "./distributionRMSD.R " + str(pRMSD) + " " + prout
    runRscript(cmdHist)



def molconvert(pfilin, pfilout= ""):
    """Convert with black background"""
    if pfilout == "":
        pfilout = pfilin[:-3] + "png"

    if path.exists(pfilout):
        return pfilout
    #cmdconvert = "molconvert \"png:w500,Q100,#00000000\" " + pfilin + " -o " + pfilout  # for transparent background
    #cmdconvert = "molconvert \"png:w500,Q100,#000000\" " + pfilin + " -o " + pfilout # black
    cmdconvert = "molconvert \"png:w500,Q100,#ffffff\" " + pfilin + " -o " + pfilout  # white
    system(cmdconvert)
    return pfilout




def babelConvertSDFtoPDB(psdf):

    cmdconvert = "babel " + psdf + " " + psdf[:-4] + ".pdb 2>/dev/null"

    if path.exists(psdf[:-4] + ".pdb") and path.getsize(psdf[:-4] + ".pdb") > 0:
        return psdf[:-4] + ".pdb"
    else:
        print cmdconvert
        system(cmdconvert)
    return psdf[:-4] + ".pdb"


def babelPDBtoMOL2(ppdbin, debug = 0):
    pfilout = ppdbin[0:-4] + ".mol2"
    if path.exists(pfilout):
        return pfilout
    else:
        cmd_convert = "babel  " + ppdbin + " " + pfilout
        #print cmd_convert
        if debug == 1:
            system(cmd_convert)
        else:
            system(cmd_convert + " 2> /dev/null")

    if not path.exists(pfilout):
        sleep(3)
        print cmd_convert
        print "SLEEP babel error"
        return babelPDBtoMOL2(ppdbin)
    else:
        return pfilout



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
    if path.exists(pcms) and path.getsize(pcms) > 0:
        return pcms

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


def multisimGDesmond(jobname, pcms, timeMDns, frameInverval, WAIT=1, HOST="gpu1"):

    prDM = path.dirname(pcms) + "/"
    pmsj = prDM + jobname + ".msj"
    pcfg = prDM + jobname + ".cfg"

    # control if out cms exisits
    pcmsout = pcms[:-4] + "-out.cms"
    if path.exists(pcmsout) and path.getsize(pcmsout) > 0:
        return "EXIST"

    # write file
    toolbox.writeFilesParamaterDesmond(pmsj, pcfg, timeMDns, frameInverval)

    # move on run folder
    prsource = getcwd()
    prrun = path.dirname(pcms)
    chdir(prrun)


    cmdDesmond = MULTISMIM + " -JOBNAME " + str(jobname) + " -maxjob 1 -cpu 1 -m " + str(pmsj) + " -c " + str(pcfg) + " " + str(pcms) + " -mode umbrella -set stage[1].set_family.md.jlaunch_opt=[\"-gpu\"] -o " + pcms[:-4] + "-out.cms -HOST " + str(HOST)

    if WAIT == 1:
        cmdDesmond = cmdDesmond + " -WAIT"
    print cmdDesmond
    system(cmdDesmond)
    chdir(prsource)

    return "DONE"

def concateneStructure(pPDBprot, pligSDF, pmaeComplex):

    if path.exists(pmaeComplex) and path.getsize(pmaeComplex) > 0:
        return pmaeComplex
    cmdStructcat = STRUCTCAT + " " + str(pPDBprot) + " " + str(pligSDF) + " -omae " + str(pmaeComplex)

    #print cmdStructcat
    system(cmdStructcat)

    if path.exists(pmaeComplex) and path.getsize(pmaeComplex) > 0:
        return pmaeComplex
    else:
        dddd
        return "ERROR"


def centerMD(ppcms, ptrj, wait = 0):

    # control existance
    pouttrj = ptrj[0:-4] + "_center_trj"
    poutcms = ptrj[0:-4] + "_center-out.cms"

    #try: removedirs(pouttrj)
    #except: pass
    #try: remove(poutcms)
    #except: pass
    #return [poutcms, pouttrj]

    #print ptrj[0:-4]
    if path.exists(pouttrj) and path.exists(poutcms):
        return [poutcms, pouttrj]
    else:
        if wait == 1:
            cmd = RUN + " -FROM desmond center.py -t " + ptrj + " " + ppcms + " " + ptrj[0:-4] + "_center"
        else:
            cmd = RUN + " -FROM desmond center.py -t " + ptrj + " " + ppcms + " " + ptrj[0:-4] + "_center&"
        print cmd
        system(cmd)
    return [poutcms, pouttrj]


def extractFrame(ppcms, ptrj, prframes, noHOH =1, step=10, MDtime=15000):

    #control numbr of frame extracted
    nbframeth = int(float(MDtime) / (step*10))
    nbframe = len(listdir(prframes))
    #print nbframeth, nbframe
    if nbframe >= nbframeth:
        #print "l.302 - cut"
        return prframes
    else:
        if noHOH == 1:
            cmd = RUN + " -FROM desmond trajectory_extract_frame.py " + str(ppcms) + " " + str(ptrj) + " -f '::" + str(step) + "' -o pdb -b " + str(prframes) + "frame 2>/dev/null&"
        else:
            cmd = RUN + " -FROM desmond trajectory_extract_frame.py " + str(ppcms) + " " + str(ptrj) + " -f '::" + str(
                step) + "' -o pdb -a 'not water' -b " + str(prframes) + "frame 2>/dev/null&"
        print "Extract Frames =>", cmd
        system(cmd)
        return prframes



def runTMalign(ppr1, ppr2, prout, debug=1):

    pathFolder.createFolder(prout)

    spdb1 = PDB.PDB(ppr1, hydrogen=1)
    spdb1.removeChain()
    print len(spdb1.latom)
    ppr1 = spdb1.writePDB(prout + path.basename(ppr1))
    print ppr1


    spdb2 = PDB.PDB(ppr2, hydrogen=1)
    spdb2.removeChain()
    ppr2 = spdb2.writePDB(prout + path.basename(ppr2))

    cmd_run = TMalign + " " + str(ppr1) + " " + str(ppr2) + " -o " + prout + "align.out -m " + prout + "matrix.out" + " > " + prout + "RMSD"
    if debug:
        print cmd_run
    system(cmd_run)

    return [prout + "align.out", prout + "align.out_all", prout + "align.out_atm",
            prout + "align.out_all_atm", prout + "matrix.out", prout + "RMSD"]


def runShaep(p_struct1, p_struct2, p_out, clean=0):

    print p_struct1
    print p_struct2
    print p_out


    if clean == 1:
        if path.exists(p_out):
            remove(p_out)
        else:
            pass
    elif path.exists(p_out):
        return p_out

    # run
    cmd = SHAEP + " --output-file " + p_out + " " + p_struct1 + " " + p_struct2 + " --noOptimization"
    print cmd
    system(cmd)
    cmd_rm = "rm " + p_out[0:-4] + "_hits.txt"

    try:
        system(cmd_rm)
    except:
        pass

    return p_out



def runscatterplotRMSD(pfilout):

    cmd = "./plotRMSD.R " + str(pfilout)
    runRscript(cmd)


def runScatterplotRMSF(pfilin):

    cmd = "./plotRMSF.R " + pfilin
    runRscript(cmd)


def scatterplotShaEP(pfilin):

    cmd = "./plotShaEP.R " + pfilin
    runRscript(cmd)


def RMSFLig(pfilin):

    cmd = "./plotRMSFlig.R " + pfilin
    runRscript(cmd)



def runRadi(ppocketatom, prout, debug=0):

    pparam = writeParameterRadi(prout, ppocketatom)
    pfilout = prout + "radi.out"

    cmd_radi = PRADI + " < " + pparam + " > " + pfilout
    if debug: print cmd_radi
    system(cmd_radi)
    system("rm " + pparam)
    return pfilout



def writeParameterRadi(prout, ppocketatom, extention="PDB"):
    """
    Write parameter file for radi
    args: -> path directory descriptors
          -> path pocket atom
          -> type file extention
    return: -> path file parameters
    """

    pfilout = prout + "parameter.radi"
    filout = open(pfilout, "w")
    filout.write(extention + "\n" + ppocketatom + "\nEPSTAB\n")
    filout.close()
    return pfilout


def runFreeSASA(ppdbin, pfilout, rsa=0):

    # asa - replace by bfactor
    cmd = "freesasa --print-as-B-value " + str(ppdbin) + " --format=pdb -o " + pfilout

    if not path.exists(pfilout):
        print cmd
        system(cmd)

    if rsa == 1:
        pfiloutrssa = pfilout[0:-4] + ".rsa"
        cmdrsa = "freesasa " + str(ppdbin) + " --format=rsa -o " + pfiloutrssa
        print cmdrsa
        system(cmdrsa)
        return [pfilout, pfiloutrssa]

    return [pfilout]



def clusterize(pdesc, pAff, typeaff, cutoff, prout):

    cmd = "./clusterDesc.R " + str(pdesc) + " " + str(pAff) + " " + str(typeaff) + " " + str(cutoff) + " ward.D2 euclidean hclust " + str(prout)
    runRscript(cmd)


def activityCliff(pdesc, paff, typeAff, Dcutoff, prout):

    cmd = "./activityCliff.R " + pdesc + " " + paff + " " + typeAff + " " + str(Dcutoff) + " ward.D2 euclidean hclust " + prout
    runRscript(cmd)




### FOR MATRIX OF DESC ###
##########################

def prepareMatrixDesc(pdesc, corcoef, maxQuantile, prout):


    cmd = "./DataPrep.R " + pdesc + " " + str(corcoef) + " " + str(maxQuantile) + " " + prout
    runRscript(cmd)

    return pdesc + "_clean.csv"

#### FOR QSAR ####
##################

def createSetFromTable(pdescglobal, ptrain, ptest, paff, prout, corcoef, maxQuantile, logAff=0, typeAff='All', nbNA=100):

    cmd = "./splitDatawithID.R %s %s %s %s %s %s %s %i %s %i"%(pdescglobal, ptrain, ptest, paff, prout, corcoef, maxQuantile, logAff, typeAff, nbNA)
    runRQSARModeling(cmd)
    dfile = {}
    if path.exists(prout + "trainSet.csv") and path.exists(prout + "testSet.csv"):
        dfile["train"] = prout + "trainSet.csv"
        dfile["test"] = prout + "testSet.csv"
    return dfile


def prepareDataset(pdesc, paff, prout, corcoef, maxQuantile, valSplit, typeAff="All", logaff=0, nbNA = 100):

    # extract train and test file
    dfile = {}
    lfile = listdir(prout)
    for filedir in lfile:
        if search("^trainSet.csv", filedir):
            dfile["train"] = prout + filedir

        elif search("^testSet.csv", filedir):
            dfile["test"] = prout + filedir

    if dfile == {}:
        cmd = "./QSARsPrep.R " + str(pdesc) + " " + str(paff) + " " + prout + " " + str(corcoef) + " " + str(
            maxQuantile) + " " + str(valSplit) + " " + str(logaff) + " " + str(typeAff) + " " + str(nbNA)
        runRQSARModeling(cmd)
        return prepareDataset(pdesc, paff, prout, corcoef, maxQuantile, valSplit, typeAff, logaff, nbNA)
    else:
        return dfile




def QSARsReg(ptrain, ptest, pcluster, prout, internalCV = 0, nbfold=10):

    cmd_QSAR = "./QSARsReg.R " + ptrain + " " + ptest + " " + pcluster + " " + prout + " " + str(nbfold) + " " + str(internalCV) + " >" + prout + "perf.txt"
    runRQSARModeling(cmd_QSAR)

    return prout + "perf.txt"


def runRQSARModeling(cmd):

    workdir = getcwd()
    chdir("/home/borrela2/development/QSARPR/source/")
    print(cmd)
    system(cmd)
    chdir(workdir)


def runRscript(cmd, out=0):

    chdir("./../Rscripts/")
    print cmd
    if out == 0:
        system(cmd)
        output = 0
    else:
        import subprocess
        output = subprocess.check_output(cmd, shell=True)
    chdir("./../py/")
    return output


def QSARsVisuData(pdescglobal, paff, prout, corcoef, maxQuantile, logAff):

    cmdVisu = "./QSARsVisuData.R " + str(pdescglobal) + " " + str(paff) + " " + str(prout) + " " + str(corcoef) + " " +\
              str(maxQuantile) + " " + str(logAff)

    runRscript(cmdVisu)


def computeRegPerf(ppred, name, prout):


    cmd = "./PerfReg.R " + ppred + " " + name + " " + prout
    runRQSARModeling(cmd)