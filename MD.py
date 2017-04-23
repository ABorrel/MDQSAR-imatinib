from os import path, listdir, makedirs

import PDB
import runExternalSoft
import toolbox


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

        nameProt = pProt.split("/")[-1][0:-4]

        dinit = {}
        for lig in llig[:30]: # add change
            namecpd = lig.split(".")[0]
            print namecpd, nameProt
            nameFolder = namecpd + "_" + nameProt
            pfolder = self.prDM + nameFolder + "/"
            if not path.exists(pfolder):
                makedirs(pfolder)

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
        step = 1
        lMDfolder = []
        while i < nbMD:
            jobname = self.lMD.keys()[i]
            #system builder
            pcms = runExternalSoft.multisimSystemBuilder(jobname, self.lMD[jobname]["complexMAE"], WAIT=1)
            self.lMD[jobname]["pcms"] = pcms
            lMDfolder.append(path.dirname(pcms) + "/")
            print len(lMDfolder)
            #GDESMON
            runExternalSoft.multisimGDesmond(jobname, pcms, self.MDtime, self.interval, WAIT=0)# add a existance criteria
            lMDfolder = toolbox.parralelLaunch(lMDfolder, self.stepWait)# control number of parralel job
            i += 1
