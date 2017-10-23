from pathFolder import listdir
from shutil import rmtree, copyfile
from numpy import mean, std

import pathFolder
import runExternalSoft
import pocketDescriptors


class MDdescriptors:

    def __init__(self, jobname, prlig, prBSs, prframe, prout):

        self.prout = pathFolder.createFolder(prout + jobname + "/")
        self.jobname = jobname
        self.prlig = prlig
        self.prBSs = prBSs
        self.prframe = prframe



    def computeLigDesc(self):

        prtemp = pathFolder.createFolder(self.prout + "ligtemp/", clean=1)

        llig = listdir(self.prlig)
        for lig in llig:
            # convert sdf
            pligsdf = runExternalSoft.babelConvertPDBtoSDF(self.prlig + lig, prtemp + lig[0:-4] + ".sdf")
        pdescframe = runExternalSoft.runPadel(prtemp, self.prout + "desc_padel", force=1)

        dout = {}
        filin = open(pdescframe, "r")
        llframeDes = filin.readlines()
        filin.close()

        ldesc = llframeDes[0].strip().split(",")

        for frameDes in llframeDes[1:]:
            i = 1  # remove name frame
            nbDesc = len(ldesc)
            ldescframe = frameDes.strip().split(",")
            while i < nbDesc:
                if not ldesc[i] in dout.keys():
                    dout[ldesc[i]] = []
                dout[ldesc[i]].append(float(ldescframe[i]))
                i += 1

        # compute desc
        for desc in dout.keys():
            av = mean(dout[desc])
            sd = std(dout[desc])

            dout[desc] = [av, sd]
            print desc, av, sd

        rmtree(prtemp)
        self.descLig = dout



    def computeBSDesc(self):



        prtemp = pathFolder.createFolder(self.prout + "BsTemp/", clean=1)

        lBSs = listdir(self.prBSs)

        for BS in lBSs:
            pBS = self.prBSs + BS, prtemp + BS
            pframe = self.prframe + .....
            copyfile(pBS)
            pocketDescriptors.pocket(ppocket, pPDB)


        for BS in lBSs:



            return



