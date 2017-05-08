from os import path, makedirs, listdir
from shutil import move, copyfile
from numpy import average

import runExternalSoft
import parseTMalign
import pathFolder
import PDB
import calculate
import toolbox
import parseShaep



class trajectoryAnalysis:

    def __init__(self, prMD, nbframe):
        self.prMD = prMD
        self.pranalyse = prMD + "analysis/"
        self.prframe = prMD + "framesMD/"
        self.nbframe = nbframe

        pathFolder.createFolder(self.pranalyse, clean=0)


    def Superimposed(self):

        self.prRMSD = self.pranalyse + "RMSDs/"
        self.matrixTranslocRot = self.prRMSD + "matrix/"
        prtemp = self.prRMSD + "temp/"

        pathFolder.createFolder(self.prRMSD)
        pathFolder.createFolder(self.matrixTranslocRot)
        pathFolder.createFolder(prtemp)

        pTMalignScore = self.prRMSD + "TMalignScore"
        if not path.exists(pTMalignScore):
            TMalignScore = open(pTMalignScore, "w")
            TMalignScore.write("Frame\tRMSD\tTMscore\n")
        else:
            TMalignScore = open(pTMalignScore, "a")

        i = 1
        imax = float(self.nbframe)
        imax = int(imax)
        while i < imax:
            print i, imax
            nframe1 = str("%05d" % (0))
            nframe2 = str("%05d" % (i))
            pframe1 = self.prframe + "frame_" + nframe1 + ".pdb"
            pframe2 = self.prframe + "frame_" + nframe2 + ".pdb"

            # control if existing
            pmatrix = self.matrixTranslocRot + str(nframe1) + "_" + str(nframe2)
            if path.exists(pmatrix) and path.getsize(pmatrix) > 0:
                i += 1
                continue
            else:
                pathFolder.createFolder(prtemp, clean=1) # clean folder temp
                print pframe1
                print pframe2

                lsuperimposed = runExternalSoft.runTMalign(pframe1, pframe2, prtemp)
                dalign = parseTMalign.parseOutputTMalign(lsuperimposed[-1])
                TMalignScore.write(str(nframe2) + "\t" + str(dalign["RMSD"]) + "\t" + str(dalign["TMscore1"]) + "\n")
                move(lsuperimposed[-2], self.matrixTranslocRot + str(nframe1) + "_" + str(nframe2))
                i += 1


        TMalignScore.close()


    def ligAnalysis(self, ligName, RMSD=1, shaepScore=0):

        self.nameLig = ligName

        if RMSD == 1:
            self.prligRMSD = self.prRMSD + "ligMD/"
            pathFolder.createFolder(self.prligRMSD)

        if shaepScore == 1:
            self.prligSheap = self.pranalyse + "ligShaep/"
            pathFolder.createFolder(self.prligSheap)
            prtempSheap = self.prligSheap + "temp/"
            pathFolder.createFolder(prtempSheap, clean=1)


        pframe0 =  self.prframe + "frame_" + str("%05d" % (0)) + ".pdb"
        cframe0 = PDB.PDB(pframe0, hydrogen=1)
        latomliginit = cframe0.get_lig(ligName)

        if RMSD == 1:
            pfiloutRMSD = self.prligRMSD + str(ligName) + "_RMSD"
            filoutRMSD = open(pfiloutRMSD, "w")
            filoutRMSD.write("Frame\tRMSD\tdistMax\n")

        if shaepScore == 1:
            pfiloutSheap = self.prligSheap + str(ligName) + "_Shaep"
            filoutSheap = open(pfiloutSheap, "w")
            filoutSheap.write("Frame\tESPscore\tshape\n")

        nb_frame = int(self.nbframe)
        nb_frame = 10
        i = 1
        dRMSFatom = {}
        while i < nb_frame:

            print "Process MD lig", i, nb_frame
            pframei = self.prframe + "frame_" + str("%05d" % (i)) + ".pdb"
            pmatrixi = self.matrixTranslocRot + str("%05d" % (0)) + "_" + str("%05d" % (i))

            cframe = PDB.PDB(pframei, hydrogen=1)
            latomligi = cframe.get_lig(ligName)
            # apply matrix transloc and rot
            matrixload = toolbox.loadMatrixTMalign(pmatrixi)
            for atomligi in latomligi:
                atomligi.applyMatrixRotTransloc(matrixload)
                if RMSD == 1:
                    if not atomligi.name in dRMSFatom.keys():
                        dRMSFatom[atomligi.name] = []
                    dRMSFatom[atomligi.name].append(atomligi)

            if shaepScore == 1:
                # wrtie frame 0
                pathFolder.cleanFolder(prtempSheap)
                pframe0Sheap = prtempSheap + "frame_" + str("%05d" % (0)) + ".pdb"
                cframe0.writePDB(pframe0Sheap, latomliginit)
                pframeisheap = prtempSheap + "frame_" + str("%05d" % (i)) + ".pdb"
                cframe.writePDB(pframeisheap, latomligi)

                runExternalSoft.runShaep(pframe0Sheap, pframeisheap, prtempSheap + "shaep-out")
                doutShaep = parseShaep.parseOutputShaep(prtempSheap + "shaep-out")
                filoutSheap.write("frame_" + str("%05d" % (i)) + "\t" + str(doutShaep["ESP_similarity"]) + "\t" + str(doutShaep["shape_similarity"]) + "\n")

            if RMSD == 1:
                lRMSD = calculate.RMSDTwoList(latomliginit, latomligi)
                filoutRMSD.write("frame_" + str("%05d" % (i)) + "\t" + str(lRMSD[0]) + "\t" + str(lRMSD[2]) + "\n")

            i += 1

        # close file
        if RMSD == 1:
            filoutRMSD.close()
        if shaepScore == 1:
            filoutSheap.close()

        # RMSF by atom for ligand
        if RMSD == 1:
            pRMSF = self.prligRMSD + self.nameLig + "_RMSF"
            filoutRMSF = open(pRMSF, "w")
            for atomname in dRMSFatom:
                filoutRMSF.write(str(atomname) + "\t")
                ldist = []
                for atomframe in dRMSFatom[atomname]:
                    for atomint in latomliginit:
                        if atomint.name == atomname:
                            ldist.append(calculate.RMSDTwoList([atomframe], [atomint])[0])
                            break
                filoutRMSF.write(str(average(ldist)) + "\n")
            filoutRMSF.close()

