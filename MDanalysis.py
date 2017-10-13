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
import superimposed


class trajectoryAnalysis:
    def __init__(self, dMD, MDtime, timeframe, stepFrame):
        self.dMD = dMD
        self.stepFrame = stepFrame
        self.MDtime = MDtime
        self.timeframe = timeframe



    def Superimpose(self, refFrame):

        prSuperimposed = self.dMD["prRMSD"] + "superimpose/"
        pathFolder.createFolder(prSuperimposed)
        self.prSuperMatrix = prSuperimposed

        prtemp = self.dMD["prRMSD"] + "temp/"
        pathFolder.createFolder(prtemp)

        pTMalignScore = self.dMD["prRMSD"] + "TMalignScore"
        if not path.exists(pTMalignScore):
            TMalignScore = open(pTMalignScore, "w")
            TMalignScore.write("Frame\tRMSD\tTMscore\n")
        else:
            TMalignScore = open(pTMalignScore, "a")

        i = self.stepFrame
        imax = float(self.MDtime) / float(self.timeframe)
        while i < imax:
            print i, imax
            nframe1 = str("%05d" % (refFrame))
            nframe2 = str("%05d" % (i))
            pframe1 = self.dMD["prframe"] + "frame_" + nframe1 + ".pdb"
            pframe2 = self.dMD["prframe"] + "frame_" + nframe2 + ".pdb"

            # control if existing
            pmatrix = prSuperimposed + str(nframe1) + "_" + str(nframe2)
            if path.exists(pmatrix) and path.getsize(pmatrix) > 0:
                i += self.stepFrame
                continue
            else:
                pathFolder.createFolder(prtemp, clean=1)  # clean folder temp

                lsuperimposed = runExternalSoft.runTMalign(pframe2, pframe1, prtemp)
                dalign = parseTMalign.parseOutputTMalign(lsuperimposed[-1])
                TMalignScore.write(str(nframe2) + "\t" + str(dalign["RMSD"]) + "\t" + str(dalign["TMscore1"]) + "\n")
                move(lsuperimposed[-2], pmatrix)
                i += self.stepFrame
        TMalignScore.close()


    def RMSDProt(self):

        prRMSDprot = self.dMD["prRMSD"] + "protein/"
        pathFolder.createFolder(prRMSDprot)
        if not "prSuperMatrix" in dir(self):
            self.Superimpose(0)

        # open reference frame
        nframeref = str("%05d" % (0))
        pframeref = self.dMD["prframe"] + "frame_" + nframeref + ".pdb"
        cprotref = PDB.PDB(PDB_input=pframeref)
        cprotref.get_atomProt()

        pfilout = prRMSDprot + "protRMSD"
        filout = open(pfilout, "w")
        filout.write("Time\tRMSDall\tRMSDC\tDmax\n0\t0\t0\t0\n")


        i = self.stepFrame
        imax = float(self.MDtime) / float(self.timeframe)
        while i < imax:
            nframe2 = str("%05d" % (i))
            pframe2 = self.dMD["prframe"] + "frame_" + nframe2 + ".pdb"

            cprot2 = PDB.PDB(PDB_input=pframe2)
            cprot2.get_atomProt()


            pmatrix = self.prSuperMatrix + str(nframeref) + "_" + str(nframe2)
            #apply matrix on frame 2
            matrixload = toolbox.loadMatrixTMalign(pmatrix)
            for atomprot2 in cprot2.latomProt:
                atomprot2.applyMatrixRotTransloc(matrixload)

            lRMSD = calculate.RMSDTwoList(cprotref.latomProt, cprot2.latomProt)
            filout.write("%s\t%s\t%s\t%s\n" % (i/100.0,lRMSD[0],lRMSD[1],lRMSD[2]))
            i += self.stepFrame

        filout.close()
        runExternalSoft.runscatterplotRMSD(pfilout)




    def ligRMSFShaEP(self, RMSF=1, ShaEPScore=1):

        prLig = self.dMD["prRMSD"] + "ligand/"
        pathFolder.createFolder(prLig)

        try:
            print self.prSuperMatrix
        except:
            print self.Superimpose(0)

        # open reference frame
        nframeref = str("%05d" % (0))
        pframeref = self.dMD["prLig"] + "LGD_" + nframeref + ".pdb"
        cligref = PDB.PDB(PDB_input=pframeref)
        cligref.get_lAtoms()


        if RMSF == 1:
            pfiloutRMSF = prLig + "ligRMSF"
            filoutRMSF = open(pfiloutRMSF, "w")
            filoutRMSF.write("Atom\tRMSF\n")


        if ShaEPScore == 1:
            pfiloutShaEP = prLig + "ligShaEP"
            filoutShaEP = open(pfiloutShaEP, "w")
            filoutShaEP.write("Time\tESPscore\tShape\n")

            prtemp = prLig + "temp/"
            pathFolder.createFolder(prtemp, clean=1)


        dRMSF = {}

        i = self.stepFrame
        imax = float(self.MDtime) / float(self.timeframe)
        while i < imax:
            print i, imax, "frame"
            nframe2 = str("%05d" % (i))
            pframe2 = self.dMD["prLig"] + "LGD_" + nframe2 + ".pdb"

            clig2= PDB.PDB(PDB_input=pframe2)
            clig2.get_lAtoms()

            pmatrix = self.prSuperMatrix + str(nframeref) + "_" + str(nframe2)
            # apply matrix on frame 2
            matrixload = toolbox.loadMatrixTMalign(pmatrix)
            for atomlig2 in clig2.latom:
                atomlig2.applyMatrixRotTransloc(matrixload)

            if RMSF == 1:
                print len(cligref.latom), len(clig2.latom)
                for atomRef in cligref.latom:
                    name = str(atomRef.name)
                    for atomframe in clig2.latom:
                        if name == atomframe.name:
                            RMSD = calculate.RMSDTwoList([atomRef], [atomframe])
                            if not name in dRMSF.keys():
                                dRMSF[name] = []
                            dRMSF[name].append(RMSD[0])
                            break

            if ShaEPScore == 1:
                pathFolder.cleanFolder(prtemp)
                ptempref = prtemp + pframeref.split("/")[-1]
                ptempframe2 = prtemp + pframe2.split("/")[-1]
                copyfile(pframeref, ptempref)
                clig2.writePDB(ptempframe2, conect=1)
                runExternalSoft.runShaep(ptempref, ptempframe2, prtemp + "shaep-out")
                doutShaep = parseShaep.parseOutputShaep(prtemp + "shaep-out")
                filoutShaEP.write(str(i/100.0) + "\t" + str(doutShaep["ESP_similarity"]) + "\t" + str(
                    doutShaep["shape_similarity"]) + "\n")
                pathFolder.cleanFolder(prtemp)

            i += self.stepFrame

        if ShaEPScore == 1:
            filoutShaEP.close()
            runExternalSoft.scatterplotShaEP(pfiloutShaEP)

        if RMSF == 1:
            for natom in dRMSF.keys():
                linew = str(natom) + "\t" + str(average(dRMSF[natom])) + "\n"
                filoutRMSF.write(linew)
            filoutRMSF.close()
            runExternalSoft.RMSFLig(pfiloutRMSF)


    def protResRMSF(self):

        prResidues = self.dMD["prRMSD"] + "residues/"
        pathFolder.createFolder(prResidues)

        try:
            print self.prSuperMatrix
        except:
            print self.Superimpose(0)

        # open reference frame
        nframeref = str("%05d" % (0))
        pframeref = self.dMD["prframe"] + "frame_" + nframeref + ".pdb"
        cprotref = PDB.PDB(PDB_input=pframeref)
        cprotref.get_byres()

        pfilout = prResidues + "resRMSD"
        filout = open(pfilout, "w")
        filout.write("NameRes\tall\tCa\tDmax\n")



        dRMSFres = {}
        i = self.stepFrame
        imax =  float(self.MDtime) / float(self.timeframe)
        while i < imax:
            print i, imax, "frame"
            nframe2 = str("%05d" % (i))
            pframe2 = self.dMD["prframe"] + "frame_" + nframe2 + ".pdb"

            cprot2 = PDB.PDB(PDB_input=pframe2)
            cprot2.get_lAtoms()

            pmatrix = self.prSuperMatrix + str(nframeref) + "_" + str(nframe2)
            # apply matrix on frame 2
            matrixload = toolbox.loadMatrixTMalign(pmatrix)
            for atomprot2 in cprot2.latom:
                atomprot2.applyMatrixRotTransloc(matrixload)
            cprot2.get_byres()

            for resname in cprot2.byres:
                res = resname.split("_")[0]
                if res in PDB.LRESSHORT:
                    numres = int(resname.split("_")[1])
                    if not numres in dRMSFres.keys():
                        dRMSFres[numres] = {}
                        dRMSFres[numres]["all"] = []
                        dRMSFres[numres]["Ca"] = []
                        dRMSFres[numres]["Dmax"] = []
                    RMSDRes = calculate.RMSDTwoList(cprotref.byres[resname], cprot2.byres[resname])
                    dRMSFres[numres]["all"].append(RMSDRes[0])
                    dRMSFres[numres]["Ca"].append(RMSDRes[1])
                    dRMSFres[numres]["Dmax"].append(RMSDRes[3])
            i += self.stepFrame

        orderednum = sorted(dRMSFres.keys())
        for num in orderednum:
            filout.write(str(num) + "\t" + str(average(dRMSFres[num]["all"])) + "\t" + str(average(dRMSFres[num]["Ca"])) + "\t" + str(average(dRMSFres[num]["Dmax"])) + "\n")
        filout.close()


        runExternalSoft.runScatterplotRMSF(pfilout)