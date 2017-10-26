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

    def __init__(self, pframe, plig, pBS, prout):

        self.pframe = pframe
        self.plig = plig
        self.pBS = pBS
        self.prout = prout


    def computeFPI(self, clean=0):


        frameID = self.plig.split("/")[-1].split("_")[1].split(".")[0]

        pfileFPI = self.prout + "FPI_" + frameID + ".csv"

        # if file exsit => load PFI
        if path.exists(pfileFPI) and clean==0:
            return pfileFPI


        fileFPI = open(pfileFPI, "w")
        dout = {}
        # header
        fileFPI.write("Ligand and pocket res\tList residues in pocket\tFPI\n")

        # define residue pocket
        cPocket = PDB.PDB(self.pBS, hydrogen=1)
        cPocket.get_byres(onlyres=1)
        lres = cPocket.getListResForFPI()

        pyplif.get_FPI(pligPDB=self.plig, ppocketPDB=self.pBS, lres=lres, filout=fileFPI)

        return pfileFPI

        #for resBS in cPocket.byres.keys():
        #    # Res as a ligand
        #    latomres = deepcopy(dres[resConsidered])
        #    PDB.changeRecoder(latomres, "HETATM")
        #    pligand = prpocket + resConsidered + ".pdb"
        #    self.CPDB.writePDB(latoms=latomres, pfilout=pligand)
        #    lplig.append(pligand)

            # Define pocket
        #    latomspocket = self.CPDB.get_BSfromLatom(latomin=latomres)
        #    ppocket = prpocket + "pocket_" + resConsidered + ".pdb"
        #    self.CPDB.writePDB(latoms=latomspocket, pfilout=ppocket)
        #    lppocket.append(ppocket)

            # format list of residues considered
        #    lresformated = PDB.convert_ListAtomtoList(latomspocket)

            # run FPI
        #    dout[resConsidered] = pyplif.get_FPI(pligPDB=pligand, ppocketPDB=ppocket, lres=lresformated, filout=fileFPI)

        #fileFPI.close()
        #self.pfileFPI = pfileFPI
        #self.FPI = dout




########################
#### FPI Analysis  #####
########################


class FPIMD:

    def __init__(self, lpFPI, jobname, prout):
        self.lpFPI = lpFPI
        self.prout = prout
        self.jobname = jobname

    def loadFPIs(self):

        dout = {}
        lresBS = []
        for pFPI in self.lpFPI:
            filin = open(pFPI, "r")
            llinesFPI = filin.readlines()
            filin.close()

            lres = llinesFPI[1].split("\t")[1].split("-")
            lFPI = llinesFPI[1].strip().split("\t")[2].split("-")
            frame = llinesFPI[1].split("\t")[0]
            dout[frame] = {}

            i = 0
            nbres = len(lFPI)
            while i < nbres:
                if lFPI[i] == "0000000":
                    i += 1
                    continue
                else:
                    if not lres[i] in lresBS:
                        lresBS.append(lres[i])
                    dout[frame][lres[i]] = lFPI[i]
                i += 1

        # write matrix FPI
        pfilout = self.prout + "FPImatrix"
        filout = open(pfilout, "w")

        filout.write("frame\t" + "\t".join(lresBS) + "\n")

        lframes = sorted(dout.keys())
        for frame in lframes:
            filout.write(str(frame))
            for resBS in lresBS:
                if not resBS in dout[frame].keys():
                    dout[frame][resBS] = "0000000"
                filout.write("\t" + str(dout[frame][resBS]))
            filout.write("\n")
        self.MDFPI = dout


    def buildTanimotoMatrix(self, mpeg=1): # need to change

        pTanimotoMatrix = self.prout + "FPITanimoto"
        TanimotoMatrix = open(pTanimotoMatrix, "w")

        lframes = sorted(self.MDFPI.keys())
        TanimotoMatrix.write("\t".join(lframes) + "\n")

        for frame in lframes:
            TanimotoMatrix.write(frame)
            for frame2 in lframes:
                bitframe1 = ""
                bitframe2 = ""
                for res in self.MDFPI[frame].keys():
                    bitframe1 = bitframe1 + self.MDFPI[frame][res]
                    bitframe2 = bitframe2 + self.MDFPI[frame2][res]
                scoreTanimoto = calculate.jaccardIndex(bitframe1, bitframe2)
                TanimotoMatrix.write("\t" + str(scoreTanimoto))
            TanimotoMatrix.write("\n")
        TanimotoMatrix.close()
        self.pTmatrix = pTanimotoMatrix



    def MDFPIbyRes(self):

        dbyres = {}
        for frame in self.MDFPI.keys():
            for res in self.MDFPI[frame].keys():
                if not res in dbyres.keys():
                    dbyres[res] = {}
                    dbyres[res]["bits"] = []
                dbyres[res]["bits"].append(self.MDFPI[frame][res])


        for res in dbyres.keys():
            dbyres[res]["TanimotoScore"] = []
            dbyres[res]["FPIAv"] = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]

            i = 0
            nFPI = len(dbyres[res]["bits"])
            while i < nFPI:
                b = 0
                while b < len(dbyres[res]["bits"][i]):
                    dbyres[res]["FPIAv"][b] += int(dbyres[res]["bits"][i][b])
                    b += 1

                j = i + 1
                while j < nFPI:
                    Tanimotores = calculate.jaccardIndex(dbyres[res]["bits"][i], dbyres[res]["bits"][j])
                    dbyres[res]["TanimotoScore"].append(Tanimotores)
                    j += 1
                i += 1

        pfiloutTanimotoRes = self.prout + "FPIbyResTanimoto"
        filouTanimotoRes = open(pfiloutTanimotoRes, "w")

        pfiloutAVFPI = self.prout + "FPIbyResMean"
        filoutAvFPI = open(pfiloutAVFPI, "w")


        filouTanimotoRes.write("MDID\t" + "\t".join(dbyres.keys()) + "\n" + self.jobname)
        filoutAvFPI.write("MDID\t" + "\t".join(dbyres.keys()) + "\n" + self.jobname)

        for res in dbyres.keys():
            filouTanimotoRes.write("\t" + str(mean(dbyres[res]["TanimotoScore"])))

            lav = []
            for bit in dbyres[res]["FPIAv"]:
                lav.append(str(bit/len(self.MDFPI.keys())))
            filoutAvFPI.write("\t" + " ".join(lav))

        filoutAvFPI.write("\n")
        filouTanimotoRes.write("\n")

        filouTanimotoRes.close()
        filoutAvFPI.close()
        self.FPIbyres = dbyres



    def DescFPI(self):

        dprofile = {}
        for frame in self.MDFPI.keys():
            profile = ""
            for res in self.MDFPI[frame].keys():
                profile = profile + self.MDFPI[frame][res]
            if not profile in dprofile.keys():
                dprofile[profile] = 0

            dprofile[profile] += 1


        linteracttype = ["Apolar", "AromFF", "AromFE", "HBA", "HBD", "ElPos", "ElNeg"]
        dinteract = {}
        for interacttype in linteracttype:
            dinteract[interacttype] = []


        for profile in dprofile.keys():
            i = 0
            dbyprofile = {}
            for interacttype in linteracttype:
                dbyprofile[interacttype] = 0
            while i < len(profile):
                dbyprofile[linteracttype[i - (7*(i/7))]] += int(profile[i])
                i +=1

            for interactiontype in dbyprofile.keys():
                dinteract[interactiontype].append(dbyprofile[interactiontype])

        print dinteract

        pfilout = self.prout + "descFPI"
        filout = open(pfilout, "w")

        # header
        filout.write("MDID\tNBprofile\tMaxIdenticProfile")
        for interacttype in linteracttype:
            filout.write("\tM" + str(interacttype) + "\tSD" + str(interacttype))
        filout.write("\n")

        filout.write(self.jobname + "\t" + str(len(dprofile.keys())) + "\t" + str(max(dprofile.values())))
        for interacttype in linteracttype:
            filout.write("\t" + str(mean(dinteract[interacttype])) + "\t" + str(std(dinteract[interacttype])))
        filout.write("\n")
        filout.close()
