from os import system, path, remove, getcwd, chdir, makedirs, listdir, removedirs
from shutil import copyfile
from re import search
from time import sleep
import subprocess

import toolbox
import PDB
import pathFolder


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


def babelConvertPDBtoSDF(ppdb, pfilout = ""):

    if pfilout == "":
        pfilout = ppdb[:-4] + ".sdf"

    cmdconvert = "/usr/bin/babel " + ppdb + " " +  pfilout + " 2>/dev/null"

    if path.exists(pfilout) and path.getsize(pfilout) > 0:
        return pfilout
    else:
        #print cmdconvert
        system(cmdconvert)
    return pfilout



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

def DescAnalysis(pdesc, paffinity, prout, valcor, PCA, corMatrix, hist, dendo, clustering):

    cmdVisu = "./visualization.R " + str(pdesc) + " " + str(paffinity) + " " + str(prout) + " " + str(valcor) + " " + \
        str(PCA) + " " + str(corMatrix) + " " + str(hist) + " " + str(dendo) + " " + str(clustering)

    print cmdVisu
    system(cmdVisu)


def MatrixMCS(pmatrix, paff, ptext, probmatrix=1):

    cmdmatrix = "./matrixMCS.R " + pmatrix + " " + paff + " " + ptext + " " + str(probmatrix)
    print cmdmatrix
    system(cmdmatrix)


def MDSMCS(pmatrix, paff):

    cmdMDS = "./MDSMCS.R " + pmatrix + " " + paff
    print cmdMDS
    system(cmdMDS)


def corPlot(pfilin, pchembl):

    cmdCor = "./corplot.R " + pfilin + " " + pchembl
    print cmdCor
    system(cmdCor)


def histAffinity(paff):

    cmdHist = "./distributionAff.R " + str(paff)
    print cmdHist
    system(cmdHist)




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
    print cmd
    system(cmd)


def runScatterplotRMSF(pfilin):

    cmd = "./plotRMSF.R " + pfilin
    print cmd
    system(cmd)


def scatterplotShaEP(pfilin):

    cmd = "./plotShaEP.R " + pfilin
    print cmd
    system(cmd)


def RMSFLig(pfilin):

    cmd = "./plotRMSFlig.R " + pfilin
    print cmd
    system(cmd)



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


#### FOR QSAR ####
##################


def prepareDataset(pdesc, paff, prout, corcoef, maxQuantile, valSplit, logaff=0):

    # extract train and test file
    dfile = {}
    lfile = listdir(prout)
    for filedir in lfile:
        if search("^train", filedir):
            dfile["train"] = prout + filedir

        elif search("^test", filedir):
            dfile["test"] = prout + filedir

    if dfile == {}:
        cmd = "./QSARsPrep.R " + str(pdesc) + " " + str(paff) + " " + prout + " " + str(corcoef) + " " + str(
            maxQuantile) + " " + str(valSplit) + " " + str(logaff)
        print cmd
        system(cmd)
        return prepareDataset(pdesc, paff, prout, corcoef, maxQuantile, valSplit)
    else:
        return dfile




def QSARsReg(ptrain, ptest, pcluster, prout, nbfold=10):

    cmd_QSAR = "./QSARsReg.R " + ptrain + " " + ptest + " " + pcluster + " " + prout + " " + str(nbfold) + " >" + prout + "perf.txt"
    print cmd_QSAR
    system(cmd_QSAR)

    return prout + "perf.txt"






