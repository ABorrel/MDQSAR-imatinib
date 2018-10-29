import pathFolder
import runExternalSoft
import toolbox

from os import path, listdir
from re import search

class sdf:


    def __init__(self, psdf, pbestpose):
        self.psdf = psdf
        if not path.exists(pbestpose):
            pathFolder.createFolder(pbestpose)
        self.pbestpose = pbestpose



    def parseSDF(self):

        lout = []
        if len(listdir(self.pbestpose)) != 0:
            lfilname = listdir(self.pbestpose)

            for filname in lfilname:
                fpose = open(self.pbestpose + filname, "r")
                rpose = fpose.read()
                fpose.close()
                dcompound = toolbox.loadPoseInDict(rpose)
                lout.append(dcompound)

        else:
            filin = open(self.psdf, "r")
            handle_read = filin.read()
            # print handle_read

            l_compound = handle_read.split("$$$$\n")[:-1]
            if len(l_compound) == 0:
                l_compound = handle_read.split("$$$$\n")

            for compound in l_compound:
                dcompound = toolbox.loadPoseInDict(compound)
                lout.append(dcompound)

        self.lc = lout
        print len(lout), "Nb pose"


    def get_dockingscore(self):
        """Keep smaller score"""

        print self.lc[0].keys()

        if not "lc" in self.__dict__:
            self.parseSDF()

        if "docking" in self.__dict__:
            return self.docking

        dscore = {}
        for compound in self.lc:
            # case where protein is included, case of XP docking
            if not "r_i_docking_score" in compound.keys():
                continue

            #print compound.keys()
            #print compound["s_m_entry_name"]
            #print compound["r_i_docking_score"]

            chemblID = compound["s_m_entry_name"].split(".")[0]
            #print chemblID

            if not chemblID in dscore.keys():
                dscore[chemblID] = {}

            if not "r_i_docking_score" in dscore[chemblID].keys():
                dscore[chemblID]["r_i_docking_score"] = float(compound["r_i_docking_score"])
                dscore[chemblID]["r_i_glide_emodel"] = float(compound["r_i_glide_emodel"])
            else:
                if float(compound["r_i_docking_score"]) < dscore[chemblID]["r_i_docking_score"]:
                    dscore[chemblID]["r_i_docking_score"] = float(compound["r_i_docking_score"])
                    dscore[chemblID]["r_i_glide_emodel"] = float(compound["r_i_glide_emodel"])

        self.docking = dscore
        return dscore



    def get_bestPose(self):

        if len(self.lc) == len(listdir(self.pbestpose)):
            print "NB best poses, previously computed:", len(self.lc)
            return



        for chemblID in self.docking:
            bestScore = self.docking[chemblID]["r_i_docking_score"]

            ipose = 0
            nbpose = len(self.lc)
            while ipose < nbpose:
                namepose = self.lc[ipose]["s_m_entry_name"]
                if chemblID == namepose.split(".")[0]:
                    if bestScore == float(self.lc[ipose]["r_i_docking_score"]):
                        pfilout = self.pbestpose + namepose + ".sdf"
                        filout = open(pfilout, "w")
                        filout.write(self.lc[ipose]["sdf"])
                        filout.close()
                        # apply a format with babel to have a proper sdf
                        runExternalSoft.babelConverttoSDF(pfilout)

                    else:
                        del self.lc[ipose]
                        nbpose = nbpose - 1
                        ipose = ipose - 1
                ipose = ipose + 1

        print "NB best poses:", len(self.lc)



    def splitPose(self, prpose):

        for pose in self.lc:
            namepose = pose["s_m_entry_name"]
            nbpose = 1
            while path.exists(prpose + namepose + "." + str(nbpose) + ".sdf"):
                nbpose = nbpose + 1

            pfilout = prpose + namepose + "." + str(nbpose) + ".sdf"
            filout = open(pfilout, "w")
            filout.write(pose["sdf"])
            filout.close()
            # apply a format with babel to have a proper sdf
            runExternalSoft.babelConverttoSDF(pfilout)

