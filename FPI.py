from os import makedirs, path
from copy import deepcopy
from numpy import mean, std

import PDB
import pyplif
import toolbox
import movie
import pathFolder
import calculate


class proteinFPI:

    def __init__(self, cPDB, prFPI):

        self.CPDB = cPDB
        self.prFPI = prFPI


    def computeFPIres(self):

        pfileFPI = self.prFPI + "FPIres.csv"
        fileFPI = open(pfileFPI, "w")
        dout = {}
        #header
        fileFPI.write("Residue Name\tList residues in pocket\tFPI\n")

        # need to define one pocket by residues -> folder of pockets
        prpocket = self.prFPI + "Pockets/"
        lppocket = []
        lplig = []
        try: makedirs(prpocket)
        except: pass

        dres = self.CPDB.get_byres()
        for resatoms in dres.keys():

            #Res as a ligand
            latomres = deepcopy(dres[resatoms])
            PDB.changeRecoder(latomres, "HETATM")
            pligand = prpocket + resatoms + ".pdb"
            self.CPDB.writePDB(latoms=latomres, pfilout=pligand)
            lplig.append(pligand)

            #Define pocket
            latomspocket = self.CPDB.get_BSfromLatom(latomin=latomres)
            ppocket = prpocket + "pocket_" + resatoms + ".pdb"
            self.CPDB.writePDB(latoms=latomspocket, pfilout=ppocket)
            lppocket.append(ppocket)

            #format list of residues considered
            lresformated = PDB.convert_ListAtomtoList(latomspocket)

            #run FPI
            dout[resatoms] = pyplif.get_FPI(pligPDB=pligand, ppocketPDB=ppocket, lres=lresformated, filout=fileFPI)


        fileFPI.close()
        self.pfileFPI = pfileFPI
        self.FPI = dout




class ligFPI:

    def __init__(self, pframe, prFPI, ligID=""):

        self.pframe = pframe
        self.CPDB = PDB.PDB(pframe)
        self.prFPI = prFPI
        self.lig = ligID

    def computeFPI(self, clean=0):

        pfileFPI = self.prFPI + "FPIlig.csv"

        # if file exsit => load PFI
        if path.exists(pfileFPI) and path.getsize(pfileFPI) > 50 and clean==0:
            dout = toolbox.loadTableFPI(pfileFPI)
            self.pfileFPI = pfileFPI
            self.FPI = dout
            return dout


        fileFPI = open(pfileFPI, "w")
        dout = {}
        # header
        fileFPI.write("Ligand and pocket res\tList residues in pocket\tFPI\n")

        # need to define one pocket by residues -> folder of pockets
        prpocket = self.prFPI + "Pockets/"
        lppocket = []
        lplig = []

        try:
            makedirs(prpocket)
        except:
            pass

        dres = self.CPDB.get_byres()
        lresconsidered = []

        for resatoms in dres.keys():
            resname = resatoms.split("_")[0]
            if self.lig != [] and self.lig != resname:
                continue
            elif resname in PDB.LRES:
                continue
            else:
                lresconsidered.append(resatoms)
                #Pocket larged -> case where we would like to analyse FPI for every res included in the BS
                #latompocketlarge = self.CPDB.get_BSfromLatom(latomin=dres[resatoms], dpocket=6, dmax=15)
                #for atompocketlarge in latompocketlarge:
                #    respocket = atompocketlarge.resName + "_" + atompocketlarge.resSeq + "_" + atompocketlarge.chainID
                #    if not respocket in lresconsidered:
                #        lresconsidered.append(respocket)

        for resConsidered in lresconsidered:
            # Res as a ligand
            latomres = deepcopy(dres[resConsidered])
            PDB.changeRecoder(latomres, "HETATM")
            pligand = prpocket + resConsidered + ".pdb"
            self.CPDB.writePDB(latoms=latomres, pfilout=pligand)
            lplig.append(pligand)

            # Define pocket
            latomspocket = self.CPDB.get_BSfromLatom(latomin=latomres)
            ppocket = prpocket + "pocket_" + resConsidered + ".pdb"
            self.CPDB.writePDB(latoms=latomspocket, pfilout=ppocket)
            lppocket.append(ppocket)

            # format list of residues considered
            lresformated = PDB.convert_ListAtomtoList(latomspocket)

            # run FPI
            dout[resConsidered] = pyplif.get_FPI(pligPDB=pligand, ppocketPDB=ppocket, lres=lresformated, filout=fileFPI)

        fileFPI.close()
        self.pfileFPI = pfileFPI
        self.FPI = dout




########################
#### FPI Analysis  #####
########################


class CompareFPIMD:

    def __init__(self, dFPI, pFPI):
        self.dFPI = dFPI
        self.pFPI = pFPI

    def MDtanimoto(self, mpeg=1): # need to change

        prout = self.pFPI + "FPItanimoto/"
        if not path.exists(prout):
            makedirs(prout)

        self.ptanimoto = prout + "FPIScore.csv"
        pfliloutPDB = prout + "protBfact_state"

        filout = open(self.ptanimoto, "w")
        filout.write("resName\tMFPI\tSDFPI\n")
        dcompare = {}
        lfilePDB = []

        i = 0
        imax = len(self.dFPI.keys())
        print imax, "l43-analysis"

        while i < (imax-1):
            j = i + 1

            print self.dFPI[self.dFPI.keys()[i]].FPI

            lresFrame = self.dFPI[self.dFPI.keys()[i]].FPI.keys()
            dcomparetemp = {}
            for resName in lresFrame:
                lbit1res = ""
                lbit2res = ""
                if not resName in dcompare.keys():
                    dcompare[resName] = []

                #retrieve list residues pocket and remove mutual res
                # case where it is not exactly the same residue
                if not resName in self.dFPI[self.dFPI.keys()[i]].FPI.keys():
                    lrespocket = self.dFPI[self.dFPI.keys()[j]].FPI[resName].keys()
                elif not resName in self.dFPI[self.dFPI.keys()[j]].FPI.keys():
                    lrespocket = self.dFPI[self.dFPI.keys()[i]].FPI[resName].keys()
                else:
                    lrespocket = self.dFPI[self.dFPI.keys()[i]].FPI[resName].keys() + self.dFPI[self.dFPI.keys()[j]].FPI[resName].keys()
                    lrespocket = list(set(lrespocket))


                for respocket in lrespocket:

                    # case frame in i
                    resbit1 = "0000000"
                    resbit2 = "0000000"

                    if resName in self.dFPI[self.dFPI.keys()[i]].FPI.keys():
                        if respocket in self.dFPI[self.dFPI.keys()[i]].FPI[resName].keys():
                            resbit1 = self.dFPI[self.dFPI.keys()[i]].FPI[resName][respocket]

                    if resName in self.dFPI[self.dFPI.keys()[j]].FPI.keys():
                        if respocket in self.dFPI[self.dFPI.keys()[j]].FPI[resName].keys():
                            resbit2 = self.dFPI[self.dFPI.keys()[j]].FPI[resName][respocket]

                    if resbit1 != "0000000" and resbit2 != "0000000":
                        lbit1res = lbit1res + resbit1
                        lbit2res = lbit2res + resbit2

                scoreJ = calculate.jaccardIndex(lbit1res, lbit2res)
                dcompare[resName].append(scoreJ)
                dcomparetemp[resName] = scoreJ

            # write PDB file
            #print self.dFPI.keys()[j]
            #print self.dFPI[self.dFPI.keys()[j]]
            #print dir(self.dFPI[self.dFPI.keys()[j]])
            #print dir(self.dFPI[self.dFPI.keys()[j]].CPDB)

            self.dFPI[self.dFPI.keys()[j]].CPDB.change_bfactor(dcomparetemp)
            self.dFPI[self.dFPI.keys()[j]].CPDB.writePDB(pfliloutPDB + str(i) + ".pdb", model=0)
            lfilePDB.append(pfliloutPDB + str(i) + ".pdb")
            i += 1


        for res in dcompare.keys():
            #lscore = dcompare[res]
            M = mean(dcompare[res])
            SD = std(dcompare[res])
            filout.write(str(res) + "\t" + str(M) + "\t" + str(SD) + "\n")
            dcompare[res] = M
        filout.close()

        if mpeg == 1:
            movie.generateMovie(lfilePDB, prout + "FPI_evolution.mpeg")


        #color by tanimoto
        self.dFPI[self.dFPI.keys()[imax-1]].CPDB.change_bfactor(dcompare)
        self.dFPI[self.dFPI.keys()[imax-1]].CPDB.writePDB(prout + "MprotBfact.pdb")


    def MDprop(self):

        prout = self.pFPI + "/FPIproba/"
        self.MFPI = prout
        pathFolder.createFolder(prout)

        self.pproba = prout + "FPIScore.csv"

        filout = open(self.pproba, "w")
        filout.write("resName\tMnbResPocket\tSDnbResPocket\tNBprofile\tMaxIdenticProfile\tMApolar\tSDApolar"
                     "\tMAromFF\tSDAromFF\tMaromFE\tSDaromFE\tMHBaccept\tSDHBaccept\tMHBdonor\tSDHBdonor"
                     "\tMElectroPos\tSDElectroPos\tMElectroNeg\tSDElectroNeg\n")

        dFPIall = {}
        dcountRes = {}

        #count in term of FPI
        dcountApolar = {}
        dcountAromFF = {}
        dcountAromEF = {}
        dcountHBa = {}
        dcountHBd = {}
        dcountElectroPos = {}
        dcountElectroNeg = {}


        i = 0
        imax = len(self.dFPI.keys())
        print imax, "NB frame (l135-analysis)"

        while i < (imax - 1):
            lresbyFrame = self.dFPI[self.dFPI.keys()[i]].FPI.keys()

            for resbyFrame in lresbyFrame:
                if not resbyFrame in dFPIall.keys():
                    dFPIall[resbyFrame] = {}
                if not resbyFrame in dcountRes.keys():
                    dcountRes[resbyFrame] = []
                    dcountApolar[resbyFrame] = []
                    dcountAromFF[resbyFrame] = []
                    dcountAromEF[resbyFrame] = []
                    dcountHBa[resbyFrame] = []
                    dcountHBd[resbyFrame] = []
                    dcountElectroPos[resbyFrame] = []
                    dcountElectroNeg[resbyFrame] = []

                lrespocket = self.dFPI[self.dFPI.keys()[i]].FPI[resbyFrame].keys()
                lrespocket = list(set(lrespocket))# order FPI


                strFPI = ""
                countResPocket = 0.0
                countApolar = 0.0
                countAromFF = 0.0
                countAromEF = 0.0
                countHBa = 0.0
                countHBd = 0.0
                countElectroPos = 0.0
                countElectroNeg = 0.0
                for respocket in lrespocket:
                    # remove FPI without interaction
                    FPIres = self.dFPI[self.dFPI.keys()[i]].FPI[resbyFrame][respocket]
                    if FPIres != "0000000":
                        strFPI = strFPI + str(respocket) + "-" + str(FPIres)
                        countResPocket += 1

                        # count by bit
                        if FPIres[0] == '1':
                            countApolar += 1
                        if FPIres[1] == '1':
                            countAromFF += 1
                        if FPIres[2] == '1':
                            countAromEF += 1
                        if FPIres[3] == '1':
                            countHBa += 1
                        if FPIres[4] == '1':
                            countHBd += 1
                        if FPIres[5] == '1':
                            countElectroPos += 1
                        if FPIres[6] == '1':
                            countElectroNeg += 1


                if not strFPI in dFPIall[resbyFrame].keys():
                    dFPIall[resbyFrame][strFPI] = 0
                dFPIall[resbyFrame][strFPI] += 1

                dcountRes[resbyFrame].append(countResPocket)
                dcountApolar[resbyFrame].append(countApolar)
                dcountAromFF[resbyFrame].append(countAromFF)
                dcountAromEF[resbyFrame].append(countAromEF)
                dcountHBa[resbyFrame].append(countHBa)
                dcountHBd[resbyFrame].append(countHBd)
                dcountElectroPos[resbyFrame].append(countElectroPos)
                dcountElectroNeg[resbyFrame].append(countElectroNeg)
            i += 1

        # write file
        for resframe in dFPIall.keys():
            filout.write(resframe)
            filout.write("\t" + str(mean(dcountRes[resframe])) + "\t" + str(std(dcountRes[resframe])))
            # stability of FPI
            ltemp = []
            for FPI in dFPIall[resframe].keys():
                ltemp.append(dFPIall[resframe][FPI])
            SumFPI = float(sum(ltemp))
            MaxFreqFPIRep = max(ltemp)
            filout.write("\t" + str(SumFPI) + "\t" + str(MaxFreqFPIRep))

            # count by bit
            filout.write("\t" + str(mean(dcountApolar[resframe])) + "\t" + str(std(dcountApolar[resframe])) + "\t" +
                         str(mean(dcountAromFF[resframe])) + "\t" + str(std(dcountAromFF[resframe])) + "\t" +
                         str(mean(dcountAromEF[resframe])) + "\t" + str(std(dcountAromEF[resframe])) + "\t" +
                         str(mean(dcountHBa[resframe])) + "\t" + str(std(dcountHBa[resframe])) + "\t" +
                         str(mean(dcountHBd[resframe])) + "\t" + str(std(dcountHBd[resframe])) + "\t" +
                         str(mean(dcountElectroPos[resframe])) + "\t" + str(std(dcountElectroPos[resframe])) + "\t" +
                         str(mean(dcountElectroNeg[resframe])) + "\t" + str(std(dcountElectroNeg[resframe])) + "\n")
        filout.close()


    def pobaFPI(self):

        if not "MFPI" in dir(self):
            prout = self.pFPI + "/FPIproba/"
            self.MFPI = prout
            pathFolder.createFolder(prout)


        dresPocket = {}
        i = 0
        imax = len(self.dFPI.keys())
        print imax, "NB frame (l135-analysis)"

        while i < imax:
            lresbyFrame = self.dFPI[self.dFPI.keys()[i]].FPI.keys()
            for resbyFrame in lresbyFrame:
                if not resbyFrame in dresPocket.keys():
                    dresPocket[resbyFrame] = []

                for resPocket in self.dFPI[self.dFPI.keys()[i]].FPI[resbyFrame].keys():
                    print
                    if not resPocket in dresPocket[resbyFrame] and self.dFPI[self.dFPI.keys()[i]].FPI[resbyFrame][resPocket] != "0000000":
                        dresPocket[resbyFrame].append(resPocket)
            i +=1


        for resConsidered in dresPocket.keys():
            dres = {}
            pfilout = self.MFPI + str(resConsidered) + ".csv"
            filout = open(pfilout, "w")

            ##########
            # header #
            ##########
            filout.write("Frame")
            for res in dresPocket[resConsidered]:
                filout.write("\t" + str(res) + "\t\t\t\t\t\t")
                dres[res] = []
            filout.write("\n")

            i=0
            imax = len(self.dFPI.keys())
            while i < imax:
                filout.write("Frame" + str(i))
                nameFrame = "frame_" + str("%05d" % (i)) + ".pdb"
                for res in dres.keys():
                    #print i
                    #print self.dFPI.keys(), "dddd"
                    #print self.dFPI.keys()[i]
                    #print resConsidered
                    #print self.dFPI[nameFrame].FPI[resConsidered]
                    #print res

                    if not res in self.dFPI[nameFrame].FPI[resConsidered].keys():
                        FPI = "0000000"
                    else:
                        FPI = self.dFPI[nameFrame].FPI[resConsidered][res]
                    dres[res].append(FPI)
                    w = "\t".join([bit for bit in FPI])
                    filout.write("\t" + w)
                filout.write("\n")
                i += 1

            filout.write("AVGFrame")
            for res in dres.keys():
                lbit = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                for FPI in dres[res]:
                    for j in range(0, 7):
                        lbit[j] += int(FPI[j])

                print len(dres[res])
                for bit in lbit:
                    filout.write("\t" + str(bit/len(dres[res])))
            filout.write("\n")
            filout.close()