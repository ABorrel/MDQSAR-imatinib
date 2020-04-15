import pathFolder
import runExternalSoft
import toolbox

from os import path, listdir
from re import search

class parseSDF:
    def __init__(self, pr_out):
        self.pr_out = pr_out

    def parseSDFFile(self, pfilin):
        self.psdf = pfilin
        lout = []
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
        print len(lout), "Nb poses"


    def parseSDFFolder(self, pr_sdf):
        """
        In case of a folder is containing several chemicals in sdf format
        """
        self.pr_sdf = pr_sdf

        lout = []
        if len(listdir(pr_sdf)) != 0:
            lfilname = listdir(pr_sdf)

            for filname in lfilname:
                fpose = open(pr_sdf + filname, "r")
                rpose = fpose.read()
                fpose.close()
                dcompound = toolbox.loadPoseInDict(rpose)
                lout.append(dcompound)
        self.lc = lout
        print len(lout), "Nb pose"


    def splitChem(self, pr_out):

        for pose in self.lc:
            namepose = pose["s_m_entry_name"]
            nbpose = 1
            while path.exists(pr_out + namepose + "." + str(nbpose) + ".sdf"):
                nbpose = nbpose + 1

            pfilout = pr_out + namepose + "." + str(nbpose) + ".sdf"
            filout = open(pfilout, "w")
            filout.write(pose["sdf"])
            filout.close()
            # apply a format with babel to have a proper sdf
            runExternalSoft.babelConverttoSDF(pfilout)

