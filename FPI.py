from os import makedirs
from copy import deepcopy

import PDB
import pyplif
import toolbox
import calculate

class proteinFPI:

    def __init__(self, cPDB, prFPI):

        self.CPDB = cPDB
        self.prFPI = prFPI


    def computeFPI(self):

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

    def __init__(self, cPDB, prFPI, ligin=""):
        """

        :param cPDB:
        :param prFPI:
        :param ligin: can be a PDB class or a lig ID
        """

        self.CPDB = cPDB
        self.prFPI = prFPI
        self.ligin = ligin
        self.BS = 0

    def computeFPI(self, BS=0):

        pfileFPI = self.prFPI + "FPIlig.csv"
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


        if type(self.ligin) == str:
            for resatoms in dres.keys():
                resname = resatoms.split("_")[0]
                if resname in PDB.LRES:
                    continue
                else:
                    lresconsidered.append(resatoms)
                    if BS == 1:
                        #Pocket larged
                        latompocketlarge = self.CPDB.get_BSfromLatom(latomin=dres[resatoms], dpocket=12, dmax=20)
                        for atompocketlarge in latompocketlarge:
                            respocket = atompocketlarge.resName + "_" + atompocketlarge.resSeq + "_" + atompocketlarge.chainID
                            if not respocket in lresconsidered:
                                lresconsidered.append(respocket)
        else:
            if BS == 1:
                latompocketlarge = self.CPDB.get_BSfromLatom(latomin=self.ligin, dpocket=12, dmax=20)
                for atompocketlarge in latompocketlarge:
                    respocket = atompocketlarge.resName + "_" + atompocketlarge.resSeq + "_" + atompocketlarge.chainID
                    if not respocket in lresconsidered:
                        lresconsidered.append(respocket)


        if BS == 1:
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
        else:
            if type(self.ligin) == str:
                latomlig = dres[self.ligin]
            else:
                latomlig = self.ligin.la

            print latomlig
            fff
            PDB.changeRecoder(latomlig, "HETATM")
            pligand = prpocket + "lig.pdb"
            self.CPDB.writePDB(latoms=latomlig, pfilout=pligand)
            lplig.append(pligand)

            # Define pocket
            latomspocket = self.CPDB.get_BSfromLatom(latomin=latomlig)
            ppocket = prpocket + "pocket.pdb"
            self.CPDB.writePDB(latoms=latomspocket, pfilout=ppocket)
            lppocket.append(ppocket)

            # format list of residues considered
            lresformated = PDB.convert_ListAtomtoList(latomspocket)

            # run FPI
            dout["lig"] = pyplif.get_FPI(pligPDB=pligand, ppocketPDB=ppocket, lres=lresformated, filout=fileFPI)


        fileFPI.close()
        self.pfileFPI = pfileFPI
        self.FPI = dout


    def compareFPI(self, dFPIin):

        if not "FPI" in dir(self):
            self.computeFPI()

        dcompare = {}
        for ligresin in self.FPI.keys():
            if not ligresin in dFPIin.keys():
                continue
            else:
                lresbit1 = ""
                lresbit2 = ""
                lrespocket = self.FPI[ligresin].key() + dFPIin[ligresin].keys()
                lrespocket = list(set(lrespocket))
                for respocket in lrespocket:
                    if not respocket in dFPIin[ligresin].keys():
                        dFPIin[ligresin][respocket] = "0000000"
                    if not respocket in self.FPI[ligresin].keys():
                        self.FPI[ligresin][respocket] = "0000000"
                    if not dFPIin[ligresin][respocket] == "0000000" and not self.FPI[ligresin][respocket] == "0000000":
                        lresbit1 = lresbit1 + self.FPI[ligresin][respocket]
                        lresbit2 = lresbit2 + dFPIin[ligresin][respocket]
            scoreJ = calculate.jaccardIndex(lresbit1, lresbit2)
            # print lresbit2
            # print lresbit1
            # print scoreJ
            dcompare[ligresin] = scoreJ
        return dcompare