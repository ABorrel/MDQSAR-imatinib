from os import path, listdir, makedirs
from pymol.preset import default_polar_contacts
from random import randint
from  multiprocessing import Process, Lock, Manager
import pprint

import PDB
import runExternalSoft
import toolbox
import MDanalysis
import pathFolder



class MD:
    def __init__(self, prDM, pranalysis, timeMD, timeframe, stepWait, nbGPU, nbCPU, stepFrame):
        self.prDM = prDM
        self.stepWait = stepWait
        self.MDtime = timeMD
        self.interval = timeframe
        self.pranalysis = pranalysis
        self.nbGPU = nbGPU
        self.nbCPU = nbCPU
        self.stepFrame = stepFrame

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

        # for MD launch
        lprframe = []
        i = 1
        for jobname in self.lMD.keys():
            print jobname, i
            if "pcmsout" in self.lMD[jobname].keys() and "prtrj" in self.lMD[jobname].keys():
                prframes = self.pranalysis + str(jobname) + "/framesMD/"
                lprframe.append(prframes)
                pathFolder.createFolder(prframes)
                nbframeth = float(self.MDtime)/(int(self.stepFrame))/10 +1
                if len(listdir(prframes)) <= int(nbframeth):# control if frame exist
                    runExternalSoft.extractFrame(self.lMD[jobname]["pcmsout"], self.lMD[jobname]["prtrj"], prframes, step=self.stepFrame, MDtime=self.MDtime)
                    lprframe = toolbox.parallelLaunch(lprframe, self.nbCPU, str(int(float(self.MDtime) / 10)))

                self.lMD[jobname]["prframe"] = prframes
            i += 1
            #print self.lMD[jobname]



    def centerFrame(self):

        lpcenter = []
        i=0
        for jobname in self.lMD.keys():
            if "pcmsout" in self.lMD[jobname].keys() and "prtrj" in self.lMD[jobname].keys():
                if len(lpcenter) >= (self.nbCPU -2):
                    loutcenterMD = runExternalSoft.centerMD(self.lMD[jobname]["pcmsout"], self.lMD[jobname]["prtrj"], wait=1)
                else:
                    loutcenterMD = runExternalSoft.centerMD(self.lMD[jobname]["pcmsout"], self.lMD[jobname]["prtrj"], wait=0)
                i += 1
                self.lMD[jobname]["pcmsout"] = loutcenterMD[0]
                self.lMD[jobname]["prtrj"] = loutcenterMD[1]
                lpcenter.append(path.dirname(self.lMD[jobname]["pcmsout"]))
                toolbox.parallelLaunch(lpcenter, self.nbCPU, "center")
                print jobname, len(lpcenter), i
            print self.lMD[jobname]





    def extractLigBSbyFrame(self, BSCutoff, namelig):

        for jobname in self.lMD.keys():
            print self.lMD[jobname]
            if "prframe" in self.lMD[jobname].keys():
                self.lMD[jobname]["prBSs"] = self.pranalysis + str(jobname) + "/BSs/"
                pathFolder.createFolder(self.lMD[jobname]["prBSs"])
                self.lMD[jobname]["prLig"] = self.pranalysis + str(jobname) + "/lig/"
                pathFolder.createFolder(self.lMD[jobname]["prLig"])

                lpframe = [self.lMD[jobname]["prframe"] + i for i in listdir(self.lMD[jobname]["prframe"])]
                nb_frame = len(listdir(self.lMD[jobname]["prframe"]))

                if len(listdir(self.lMD[jobname]["prLig"])) >= nb_frame and len(listdir(self.lMD[jobname]["prBSs"])) >= nb_frame:
                    continue
                else:

                    for pframe in lpframe:
                        cPDB = PDB.PDB(pframe)
                        latomlig = cPDB.get_lig(namelig)

                        cPDB.get_BSfromlig(dpocket=BSCutoff)
                        # add step of rename atom
                        pLGD = self.lMD[jobname]["prLig"] + "LGD_" + pframe.split("_")[-1]
                        pBS = self.lMD[jobname]["prBSs"] + "BS_" + pframe.split("_")[-1]

                        cPDB.writePDB(pLGD, latomlig, conect=1)
                        cPDB.writePDB(pBS, cPDB.pocketsRES["UNK_900_A"])# default in schrodinger



    def analyseRMSD(self):

        for jobname in self.lMD.keys():
            prRMSD = self.pranalysis + str(jobname) + "/RMSDs/"
            self.lMD[jobname]["prRMSD"] = prRMSD

            cMDanalysis = MDanalysis.trajectoryAnalysis(self.lMD[jobname], self.MDtime, self.interval, self.stepFrame)
            cMDanalysis.Superimpose(0)
            cMDanalysis.RMSDProt()
            cMDanalysis.protResRMSF()
            cMDanalysis.ligRMSFShaEP()



