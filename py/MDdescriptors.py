from pathFolder import listdir
from shutil import rmtree, copyfile
from numpy import mean, std
from os import path

import pathFolder
import runExternalSoft
import pocketDescriptors
import FPI
import toolbox
import descriptors3D


class MDdescriptors:

    def __init__(self, jobname, prlig, prBSs, prframe, prout):

        self.prout = prout
        self.jobname = jobname
        self.prlig = prlig
        self.prBSs = prBSs
        self.prframe = prframe



    def computeLigDesc(self):

        pfilout = self.prout + "descLig"
        #uncomment for run
        #if path.exists(pfilout) and path.getsize(pfilout)>200:
        #    return

        # folder temp to include sdf
        prtemp = pathFolder.createFolder(self.prout + "ligtemp/", clean=0)

        #run new script for 3D descriptors
        llig = listdir(self.prlig)
        llig = sorted(llig)# order lig

        pdescByFrames = self.prout + "Ligbyframe"
        fdescByFrames = open(pdescByFrames, "w")

        ldesc3D = descriptors3D.l3D
        fdescByFrames.write("Frame" + "\t" + "\t".join(ldesc3D) + "\n")

        for lig in llig:
            pligsdf = runExternalSoft.babelConverttoSDF(self.prlig + lig, prtemp + lig[0:-4] + ".sdf")
            ddesc = descriptors3D.get3Ddesc(pligsdf)
            fdescByFrames.write(lig[0:-4])
            for desc3D in ldesc3D:
                fdescByFrames.write("\t" + str(ddesc[desc3D]))
            fdescByFrames.write("\n")
        fdescByFrames.close()

        dout = {}
        filin = open(pdescByFrames, "r")
        llframeDes = filin.readlines()
        filin.close()

        ldesc = llframeDes[0].strip().split("\t")

        for frameDes in llframeDes[1:]:
            i = 1  # remove name frame
            nbDesc = len(ldesc)
            ldescframe = frameDes.strip().split("\t")
            while i < nbDesc:
                if not ldesc[i] in dout.keys():
                    dout[ldesc[i]] = []
                #print "////", (ldescframe[i])
                try:dout[ldesc[i]].append(float(ldescframe[i])) # remove frame where desc is not computed
                except:pass
                i += 1

        # compute desc
        for desc in dout.keys():
            av = mean(dout[desc])
            sd = std(dout[desc])

            dout[desc] = [av, sd]
            print desc, av, sd

        try:rmtree(prtemp)
        except:pass
        self.descLig = dout
        self.writeMDdesc(self.descLig, pfilout, self.jobname)


    def computeBSDesc(self):


        pfilout = self.prout + "descBS"
        if path.exists(pfilout):
            return

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


        pfilout = self.prout + "descFPI"
        if path.exists(pfilout):
            return

        prtemp = pathFolder.createFolder(self.prout + "FPITemp/", clean=1)

        llig = listdir(self.prlig)
        llig = sorted(llig)


        lpFPI = []
        # compute FPI for lig
        for lig in llig:
            frameID = lig.split("_")[1].split(".")[0]

            plig = prtemp + "LGD_" + frameID + ".pdb"
            pframe = prtemp + "frame_" + frameID + '.pdb'
            pBS = prtemp + "BS_" + frameID + ".pdb"

            copyfile(self.prlig + lig, plig)
            copyfile(self.prframe + "frame_" + frameID + ".pdb", pframe)
            copyfile(self.prBSs + "BS_" + frameID + ".pdb", pBS)

            CFPI = FPI.ligFPI(pframe, plig, pBS, prtemp)
            pFPI = CFPI.computeFPI()
            lpFPI.append(pFPI)


        cMDPFI = FPI.FPIMD(lpFPI, self.jobname, self.prout)
        cMDPFI.loadFPIs()

        # matrix tanimoto
        cMDPFI.buildTanimotoMatrix()
        cMDPFI.MDFPIbyRes()

        # desc
        cMDPFI.DescFPI()
        rmtree(prtemp)

        self.descFPI = cMDPFI




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




