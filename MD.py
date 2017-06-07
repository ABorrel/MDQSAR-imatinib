from os import path, listdir, makedirs

import PDB
import runExternalSoft
import toolbox
import MDanalysis
import pathFolder


class MD:
    def __init__(self, prDM, timeMD, timeframe, stepWait):
        self.prDM = prDM
        self.stepWait = stepWait
        self.MDtime = timeMD
        self.interval = timeframe

    def initialisation(self, prLigand, pProt):
        """Prot need to be prep before MD !!!!"""

        print "Protein need to be prep !!"

        # create folder for DM
        llig = listdir(prLigand)
        llig = ["CHEMBL207986.1.pdb", "CHEMBL384304.1.pdb"] # !!!!!!!!!!!

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


    def runMultipleMD(self, runMD = 1):


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
            #GDESMON
            if runMD == 1:
                runExternalSoft.multisimGDesmond(jobname, pcms, self.MDtime, self.interval, WAIT=0)# add a existance criteria
                self.lMD[jobname]["pcms-out"] = pcms[0:-4] + "-out.cms"
                self.lMD[jobname]["trj"] = path.dirname(pcms) + "/" + jobname + "_trj"
                lMDfolder = toolbox.parralelLaunch(lMDfolder, self.stepWait)# control number of parralel job
            else:
                if path.exists(pcms[0:-4] + "-out.cms"):
                    self.lMD[jobname]["pcms-out"] = pcms[0:-4] + "-out.cms"
                    self.lMD[jobname]["trj"] = path.dirname(pcms) + "/" + jobname + "_trj"
            i += 1




    def extractFrame(self):
        """Extract frame and wrap water"""

        # launch DM to build result structure
        self.runMultipleMD(runMD=0)

        # for MD launch
        for jobname in self.lMD.keys():
            print jobname
            if "pcms-out" in self.lMD[jobname].keys():
                prframes = path.dirname(self.lMD[jobname]["pcms-out"]) + "/framesMD/"
                pathFolder.createFolder(prframes)
                if len(listdir(prframes)) == 0:
                    centerMD = runExternalSoft.centerMD(self.lMD[jobname]["pcms-out"], self.lMD[jobname]["trj"])
                    self.lMD[jobname]["pcms-out"] = centerMD + "-out.cms"
                    self.lMD[jobname]["trj"] = centerMD + "_trj"
                    runExternalSoft.extractFrame(self.lMD[jobname]["pcms-out"], self.lMD[jobname]["trj"], prframes)
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







