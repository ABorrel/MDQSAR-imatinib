from pathFolder import listdir
from shutil import rmtree, copyfile
from numpy import mean, std

import pathFolder
import runExternalSoft
import pocketDescriptors
import FPI


class MDdescriptors:

    def __init__(self, jobname, prlig, prBSs, prframe, nameLig, prout):

        self.prout = pathFolder.createFolder(prout + jobname + "/")
        self.jobname = jobname
        self.prlig = prlig
        self.prBSs = prBSs
        self.prframe = prframe
        self.nameLig = nameLig



    def computeLigDesc(self):

        prtemp = pathFolder.createFolder(self.prout + "ligtemp/", clean=1)

        llig = listdir(self.prlig)
        llig = sorted(llig)
        for lig in llig:
            # convert sdf
            pligsdf = runExternalSoft.babelConvertPDBtoSDF(self.prlig + lig, prtemp + lig[0:-4] + ".sdf")
        pdescframe = runExternalSoft.runPadel(prtemp, self.prout + "Ligbyframe", force=1)

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
                print (ldescframe[i])
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
        self.writeMDdesc(self.descLig, self.prout + "descLig", self.jobname)


    def computeBSDesc(self):

        prtemp = pathFolder.createFolder(self.prout + "BsTemp/", clean=1)
        pfilebyframe = self.prout + "BSbyframe"
        filoutbyframe = open(pfilebyframe, "w")


        lBSs = listdir(self.prBSs)
        lBSs = sorted(lBSs)
        dout = {}
        header = 0
        for BS in lBSs:
            frameID = BS.split("_")[1].split(".")[0]

            # copy in folder
            pBS = prtemp + BS
            pframe = prtemp + "frame_" + frameID + ".pdb"

            copyfile(self.prBSs + BS, pBS)
            copyfile(self.prframe + "frame_" + frameID + ".pdb", pframe)


            cBS = pocketDescriptors.pocket(pBS, pframe)
            cBS.get_alldescs()

            if header == 0:
                ldesc = cBS.ldescs
                filoutbyframe.write("frame\t" + "\t".join(ldesc) + "\n")
                header = 1
            filoutbyframe.write(str(frameID))

            for desc in ldesc:
                if desc in cBS.compo.keys():
                    val = cBS.compo[desc]
                elif desc in cBS.energy.keys():
                    val = cBS.energy[desc]
                else:
                    val = cBS.geometry[desc]

                filoutbyframe.write("\t" + str(val))
            filoutbyframe.write("\n")



            for descenergy in cBS.energy.keys():
                if not descenergy in dout.keys():
                    dout[descenergy] = []
                dout[descenergy].append(float(cBS.energy[descenergy]))


            for desgeo in cBS.geometry.keys():
                if not desgeo in dout.keys():
                    dout[desgeo] = []
                dout[desgeo].append(float(cBS.geometry[desgeo]))


            for desccompo in cBS.compo.keys():
                if not desccompo in dout.keys():
                    dout[desccompo] = []
                dout[desccompo].append(float(cBS.compo[desccompo]))

        filoutbyframe.close()


        for descBS in dout.keys():
            av = mean(dout[descBS])
            sd = std(dout[descBS])
            dout[descBS] = [av, sd]

        rmtree(prtemp)
        self.descBS = dout

        self.writeMDdesc(self.descBS, self.prout + "descBS", self.jobname)




    def computeFPI(self):

        prtemp = pathFolder.createFolder(self.prout + "FPITemp/", clean=1)

        llig = listdir(self.prlig)
        llig = sorted(llig)

        for lig in llig:
            frameID = lig.split("_")[1].split(".")[0]

            plig = prtemp + "LGD_" + frameID + ".pdb"
            pframe = prtemp + "frame_" + frameID + '.pdb'

            copyfile(self.prlig + lig, plig)
            copyfile(self.prframe + "frame_" + frameID + ".pdb", pframe)


            CFPI = FPI.ligFPI(pframe, prtemp, ligID=self.nameLig)
            CFPI.computeFPI()
            dddd

            dcFpI[frameName] = CFPI

            FPIMD = FPI.CompareFPIMD(dcFpI, prtempMDFPI)
            FPIMD.MDprop()
            FPIMD.pobaFPI()
            # FPIMD.MDtanimoto() # useless if only ligand is considered

            return dout

        return



    def writeMDdesc(self, ddesc, pfilout, nameMD):

        filout = open(pfilout, "w")
        filout.write("MDID")

        ldesc = ddesc.keys()
        for desc in ldesc:
            filout.write("\tM_" + str(desc) + "\tSD_" + str(desc))
        filout.write("\n")

        filout.write(str(nameMD))
        for desc in ldesc:
            filout.write("\t" + str(ddesc[desc][0]) + "\t" + str(ddesc[desc][1]))
        filout.close()




