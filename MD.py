from os import path, listdir, makedirs
from random import randint

import PDB
import runExternalSoft
import toolbox
import MDanalysis
import pathFolder


class MD:
    def __init__(self, prDM, pranalysis, timeMD, timeframe, stepWait, nbGPU, nbCPU):
        self.prDM = prDM
        self.stepWait = stepWait
        self.MDtime = timeMD
        self.interval = timeframe
        self.pranalysis = pranalysis
        self.nbGPU = nbGPU
        self.nbCPU = nbCPU

    def initialisation(self, prLigand, pProt):
        """Prot need to be prep before MD !!!!"""

        print "Protein need to be prep !!"

        # create folder for DM
        llig = listdir(prLigand)

        nameProt = pProt.split("/")[-1][0:-4]

        dinit = {}
        for lig in llig: # add change
            namecpd = lig.split(".")[0]
            print namecpd, nameProt
            nameFolder = namecpd + "_" + nameProt
            pfolder = self.prDM + nameFolder + "/"

            pathFolder.createFolder(pfolder)

            #create PDB file as input
            pMAE = pfolder + "complexInit.mae"
            pligandSDF = prLigand + lig
            pprotPDB = pProt

            pMAE = runExternalSoft.concateneStructure(pprotPDB, pligandSDF, pMAE)
            if pMAE == "ERROR":
                continue
            else:
                dinit[nameFolder] = {}
                dinit[nameFolder]["pfolder"] = pfolder
                dinit[nameFolder]["complexMAE"] = pMAE# complex in input of MD
                dinit[nameFolder]["ligandSDF"] = prLigand + lig
                dinit[nameFolder]["protPrep"] = pProt

        self.lMD = dinit


    def runMultipleMD(self):


        if not "lMD" in dir(self):
            print "ERROR initialisation"
            return "ERROR - MD initialisation"

        nbMD = len(self.lMD.keys())
        i = 0
        lMDfolder = []
        while i < nbMD:
            jobname = self.lMD.keys()[i]
            #system builder
            pcms = runExternalSoft.multisimSystemBuilder(jobname, self.lMD[jobname]["complexMAE"], WAIT=1)
            self.lMD[jobname]["pcms"] = pcms
            lMDfolder.append(path.dirname(pcms) + "/")
            print len(lMDfolder)
            #GDESMOND
            if self.nbGPU != 0:
                HOSTGPU = "gpu" + str(randint(1, self.nbGPU))
                runExternalSoft.multisimGDesmond(jobname, pcms, self.MDtime, self.interval, WAIT=0, HOST=HOSTGPU)# add a existance criteria
            else:
                runExternalSoft.multisimGDesmond(jobname, pcms, self.MDtime, self.interval, WAIT=0)# add a existance criteria
            lMDfolder = toolbox.parallelLaunch(lMDfolder, self.nbGPU, "-out\.cms", self.stepWait)# control number of parralel job
            i += 1

        # recup frame
        for jobname in self.lMD.keys():
            #control if out exist
            pcmsout = self.lMD[jobname]["pfolder"] + jobname + "-out.cms"
            prtrj = self.lMD[jobname]["pfolder"] + jobname + "_trj"
            if path.exists(pcmsout) and path.getsize(pcmsout) > 100:
                self.lMD[jobname]["pcmsout"] = pcmsout

            if path.exists(prtrj):
                self.lMD[jobname]["prtrj"] = prtrj


    def extractFrame(self):
        """Extract frame and wrap water"""

        print self.lMD["CHEMBL3617738_2hyy_MD"].keys()
        print self.lMD["CHEMBL3617738_2hyy_MD"]["ligandSDF"]
        print self.lMD["CHEMBL3617738_2hyy_MD"]["pcmsout"]

        nbjreceiptboob = len(self.lMD.keys())

        # for MD launch
        for jobname in self.lMD.keys():
            print jobname
            if "pcmsout" in self.lMD[jobname].keys() and "prtrj" in self.lMD[jobname].keys():
                prframes = self.pranalysis + str(jobname) + "/framesMD/"
                pathFolder.createFolder(prframes)
                if len(listdir(prframes)) == 0:# control if frame exist
                    centerMD = runExternalSoft.centerMD(self.lMD[jobname]["pcmsout"], self.lMD[jobname]["prtrj"])
                    self.lMD[jobname]["pcmsout"] = centerMD + "-out.cms"
                    self.lMD[jobname]["prtrj"] = centerMD + "_trj"
                    runExternalSoft.extractFrame(self.lMD[jobname]["pcmsout"], self.lMD[jobname]["prtrj"], prframes)
                self.lMD[jobname]["prframe"] = prframes



    def analyseAllMD(self, RMSD=1, ligAnalysis=1, nameLig="UNK"):

        lanalysis = []
        for MDID in self.lMD.keys():
            cMDanalysis = MDanalysis.trajectoryAnalysis(self.lMD[MDID]["pfolder"], float(self.MDtime)/float(self.interval))
            cMDanalysis.Superimposed()
            if RMSD == 1:
                cMDanalysis.protResRMSF()
            if ligAnalysis == 1:
                cMDanalysis.ligAnalysis(nameLig)
            lanalysis.append(cMDanalysis)







