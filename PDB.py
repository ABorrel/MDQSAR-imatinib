'''
Created on Feb 10, 2017

@author: Alexandre Borrel
'''

from os import path, makedirs
from re import search
from math import sqrt
from copy import deepcopy
from urllib import urlretrieve
from shutil import move
import random

import runExternalSoft
import toolbox

PRPDBDATABSE = "/home/PDB/"

R = random.random() # corse random
LRES = ["ALA", "ILE", "LEU", "VAL", "MET", "CYS", "PHE", "TRP", "TYR", "HIS", "THR", "SER", "ASN", "GLN", "ASP", "GLU",
        "ARG", "LYS", "PRO", "GLY"]
#LTYPE = ["Oox", "Oh", "Oph",  "Oc", "Ow", "Nam", "Nim", "Ngu", "NaI", "Car", "Su", "Xot"]
MAXCONNECT = 2.0
MAXDISTBYRES = 12
#LRESCONSIDERD = ["HIS", "ASP", "GLU", "ARG", "LYS"] # consider all atom




class PDB:
    '''
    Load a PDB file
    '''

    def __init__(self, PDB_input, hydrogen=0):
        '''
        Build PDB object, used path
        - parse PDB with
        - option of parsing, to reduce the load (complete, coord, header)
        Note: dynamic update of PDB database
        '''


        # format input
        if path.exists(PDB_input):
            self.PDB = PDB_input[-8:-4].lower
            self.filePDB = PDB_input
        elif len(PDB_input) == 4:
            self.PDB = PDB_input
            self.filePDB = PRPDBDATABSE + str(self.PDB.lower()) + ".pdb"
        else:
            print "ERROR -> PDB input l.48 in PDB"
        #print self.filePDB
        if PDB_input[-3:] == "sdf":
            pfilePDB = runExternalSoft.babelConvertSDFtoPDB(PDB_input)
            self.filePDB = pfilePDB

        if not path.exists(self.filePDB):# case of PDB do not exisit in the folder
            if self.download_PDB() == 0:
                self.llines = []
            else:
                filePDB = open(self.filePDB, "r")
                self.llines = filePDB.readlines()
                filePDB.close()
        else:
            filePDB = open(self.filePDB, "r")
            self.llines = filePDB.readlines()
            filePDB.close()
        # keep hydrogen
        self.hydrogen = hydrogen


    def get_BSfromlig(self, lligatoms=[], dpocket = 4.5, dmax = 12):

        self.get_byres()
        self.lgd = {}
        # extract lig
        if lligatoms == []:
            for resName in self.byres.keys():
                if not resName.split("_")[0] in LRES:
                   if not resName in self.lgd.keys():
                       self.lgd[resName] = self.byres[resName]

        # case where ligand is extract from a another protein
        else:
            self.lgd["XXX"] = lligatoms
        self.pockets = {}
        self.pocketsRES = {}
        # define pocket
        for lgdname in self.lgd.keys():
            self.pockets[lgdname] = []
            self.pocketsRES[lgdname] = []

            for atomlgd in self.lgd[lgdname]:
                for res in self.byres.keys():
                    if res == lgdname:
                        continue
                    flagres = 0
                    for atomres in self.byres[res]:
                        disttest = atomres.euclidiendist(atomlgd)
                        if disttest <= dpocket:
                            if not atomres in self.pockets[lgdname]:
                                self.pockets[lgdname].append(atomres)
                                flagres = 1
                            else:
                                continue
                        elif disttest >= dmax:
                            break

                        # case where pocket is define based on residues
                        if flagres == 1:
                            for atomres in self.byres[res]:
                                if not atomres in self.pocketsRES[lgdname]:
                                    self.pocketsRES[lgdname].append(atomres)
                            flagres=0


    def get_BSfromLatom(self, latomin, dpocket=4.5, dmax=12, dmin=2, delconnected=1):

        lout = []
        lresconnected = []
        lresextracted = []
        if not "byres" in self.__dict__:
            self.get_byres()

        for atomin in latomin:
            for res in self.byres.keys():
                for atomres in self.byres[res]:
                    disttest = atomres.euclidiendist(atomin)
                    if disttest < dmin and delconnected == 1:
                        if not res in lresconnected:
                            lresconnected.append(res)
                    elif disttest <= dpocket:
                        #print disttest, "dist"
                        if not res in lresextracted:
                            lresextracted.append(res)
                    elif disttest >= dmax:
                        break

        for resextracted in lresextracted:
            if not resextracted in lresconnected:
                for atom in self.byres[resextracted]:
                    lout.append(atom)
        return lout


    def get_chains(self,lchainID):
        l_out = []
        for linePDB in self.llines:
            if search("^HEADER", linePDB) or search("^TITLE", linePDB) or search("^REMARK", linePDB):
                l_out.append(linePDB)
            elif search("^ATOM", linePDB) or search("^HETATM", linePDB):
                atom = Atom(linePDB)
                if atom.chainID in lchainID:
                    l_out.append(linePDB)
        return l_out


    def download_PDB(self, debug=1):
        urlSeq = ("http://www.pdb.org/pdb/files/%s.pdb" % self.PDB)

        try:
            ppdb = urlretrieve(urlSeq)
            if debug: print ppdb
            move(ppdb[0], PRPDBDATABSE + self.PDB + ".pdb")
            print str(self.PDB) + "-> done"
            return 1
        except:
            print str(self.PDB) + "-> ERROR DOWNLOAD PDB file (l.64)"
            return 0

    def get_resolution(self):

        try:
            return self.rX
        except:
            pass
        nb_line = len(self.llines)
        resolution = -1
        i = 0
        while i < nb_line:
            if resolution != -1:
                try:
                    self.rX = float(resolution)
                except:
                    self.rX = "NA"
                return self.rX
            if search("^REMARK   2 RESOLUTION", self.llines[i]):
                try:
                    resolution = self.llines[i].split("RESOLUTION.")[-1].split("ANGSTROM")[0].replace(" ", "")
                except:
                    pass
            i += 1

        self.rX = "NA"
        return "NA"

    def get_Rfree(self):

        try:
            return self.Rfree
        except:
            pass

        nb_line = len(self.llines)
        Rfree = -1
        i = 0
        while i < nb_line:
            if Rfree != -1:
                try:
                    self.Rfree = float(Rfree)
                except:
                    self.Rfree = "NA"
                return self.Rfree
            if search("REMARK   3   R VALUE", self.llines[i]):
                try:
                    Rfree = self.llines[i].strip().split(":")[-1].replace(" ", "")
                except:
                    pass
            i += 1

        self.Rfree = "NA"
        return "NA"

    def get_method(self):

        try:
            return self.expmeth
        except:
            pass

        nb_line = len(self.llines)
        expmeth = -1
        i = 0
        while i < nb_line:
            if expmeth != -1:
                self.expmeth = expmeth
                return expmeth
            if search("EXPDTA", self.llines[i]):
                try:
                    expmeth = self.llines[i].split("EXPDTA")[-1].replace(" ", "").strip()
                except:
                    pass
            i += 1

        self.expmeth = "NA"
        return "NA"

    def get_lAtoms(self):

        try:
            return self.latom
        except:
            pass

        l_out = []
        for linePDB in self.llines:
            # only cordinate atom considered
            if search("^ATOM", linePDB) or search("^HETATM", linePDB):
                atomadd = Atom(linePDB)
                if atomadd.element == "H" and self.hydrogen != 0:
                    l_out.append(atomadd)
                elif atomadd.element != "H":
                    l_out.append (atomadd)
        self.latom = l_out
        return l_out

    def get_byres(self):

        try:
            return self.byres
        except:
            pass

        try:
            self.latom
        except:
            self.get_lAtoms()

        d_res = {}
        for atom in self.latom:
            if atom.chainID == " ":
                k_in = str(atom.resName) + "_" + str(atom.resSeq) + "_0"
            else:
                k_in = str(atom.resName) + "_" + str(atom.resSeq) + "_" + str(atom.chainID)
            if not k_in in d_res.keys():
                d_res[k_in] = []
            d_res[k_in].append(atom)
        self.byres = d_res
        return d_res

    def get_llig(self):

        try:
            return self.llig
        except:
            pass

        l_out = []
        for atom in self.get_lAtoms():
            if not atom.resName in LRES and atom.recoder == "HETATM":
                if not atom.resName in l_out and atom.resName != "HOH":
                    l_out.append(atom.resName)

        self.llig = l_out
        return l_out

    def get_header(self):

        try:
            return self.header
        except:
            pass
        header = self.llines[0][6:].lower().strip()
        self.header = header
        return header

    def change_bfactor(self, dp_bfactor):
        """

        :param d_bfactor: dictionary with residue and atom serial bfactor (*100 on score)
        :return: NONE
        """
        if path.exists(dp_bfactor):
            d_bfactor = {}
            filin = open(dp_bfactor, "r")
            llines = filin.readlines()
            filin.close()

            for lineFPI in llines[1:]:
                lelem = lineFPI.strip().split("\t")
                print lelem
                d_bfactor[lelem[0]] = float(lelem[1])
        else:
            d_bfactor = dp_bfactor


        self.get_byres()
        for res in self.byres.keys():
            i = 0
            while i < len(self.byres[res]):
                atom = self.byres[res][i]
                k_in = atom.resName + "_" + str(atom.resSeq) + "_" + str(atom.chainID)
                if k_in in d_bfactor:
                    if d_bfactor[k_in] == 1.00:
                        atom.Bfact = 0.00
                    else :
                        try: atom.Bfact = abs(100 - float(d_bfactor[k_in]*99.0))
                        except: atom.Bfact = 0.00
                else:
                    atom.Bfact = 00.0
                i = i + 1



    def addLigand(self, plig):


        print plig
        if not "latom" in dir(self):
            self.get_lAtoms()

        clig = PDB(plig, hydrogen=1)
        latomlig = clig.get_lAtoms()
        print len(latomlig)

        print len(self.latom)
        self.latom = self.latom + latomlig
        print len(self.latom)

    def removeChain(self):

        if not "latom" in dir(self):
            self.get_lAtoms()

        for atom in self.latom:
            atom.chainID = ""

    def get_lig(self, nameLig):

        if not "byres" in self.__dict__:
            self.get_byres()

        for res in self.byres.keys():
            if search(nameLig, res):
                # !! only first out
                return self.byres[res]



    def writePDB (self, pfilout, latoms=""):
        """
        Need add header
        :param pfilout: path of file to write
        :return:
        """

        filout = open(pfilout, "w")

        if latoms == "":
            if not "latom" in dir(self):
                latoms = self.get_lAtoms()
            else:
                latoms = self.latom

        for atom in latoms:
            filout.write("%-6s%5s %4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n" %(atom.recoder, atom.serial, atom.name, atom.char, atom.resName, atom.chainID, atom.resSeq, atom.iCode, atom.x, atom.y, atom.z, atom.occupancy, atom.Bfact, atom.element, atom.charge))
        filout.close()

        return pfilout






class Atom:
    def __init__(self, linePDB):
        if linePDB != "":
            self.recoder = linePDB[0:6].replace(" ", "")
            self.serial = linePDB[6:11].replace(" ", "")
            self.name = linePDB[12:16].replace(" ", "")
            self.char = linePDB[16]
            self.resName = linePDB[17:20].replace(" ", "")
            self.chainID = linePDB[21]
            # format in A
            if self.chainID == " ":
                self.chainID = "A"
            self.resSeq = linePDB[22:26].replace(" ", "")
            self.iCode = linePDB[26]
            self.x = float(linePDB[30:38].replace(" ", ""))
            self.y = float(linePDB[38:46].replace(" ", ""))
            self.z = float(linePDB[46:54].replace(" ", ""))
            self.occupancy = float(linePDB[54:60].replace(" ", ""))
            self.Bfact = float(linePDB[60:66].replace(" ", ""))
            self.element = linePDB[76:78].replace(" ", "")
            # case of model when element is not present
            if self.element == "":
                self.element = self.name[0]
            self.charge = linePDB[78:80].replace(" ", "")
        else:
            self.name = ""
            self.serial = ""
            self.char = ""
            self.resName = ""
            self.chainID = ""
            self.resSeq = ""
            self.iCode = ""
            self.x = ""
            self.y = ""
            self.z = ""
            self.element = ""
            self.charge = ""
            self.occupancy = ""
            self.Bfact = ""
            self.recoder = ""

    def get_connectmatrice(self, latom):

        try:
            self.connect
        except:
            pass

        l_connect = [self.serial]

        for atom in latom:
            d = self.euclidiendist(atom)
            if d <= MAXCONNECT and not atom.serial in l_connect:
                l_connect.append(atom.serial)

        self.connect = l_connect

    def euclidiendist(self, atom):

        x1 = float(self.x)
        x2 = float(atom.x)
        xd = x2 - x1

        y1 = float(self.y)
        y2 = float(atom.y)
        yd = y2 - y1

        z1 = float(self.z)
        z2 = float(atom.z)
        zd = z2 - z1

        return sqrt(xd * xd + yd * yd + zd * zd)

    def add_dist(self, dist):
        self.dist = dist


    def applyMatrixRotTransloc(self, matrixin):

        if not type(matrixin) == dict:
            matrix_transloc = toolbox.loadMatrixTMalign(matrixin)
        else:
            matrix_transloc = matrixin

        atomx = matrix_transloc["t1"] + matrix_transloc["u11"] * float(self.x) + matrix_transloc["u12"] * float(
                self.y) + matrix_transloc["u13"] * float(self.z)
        atomy = matrix_transloc["t2"] + matrix_transloc["u21"] * float(self.x) + matrix_transloc["u22"] * float(self.y) + \
                    matrix_transloc["u23"] * float(self.z)
        atomz = matrix_transloc["t3"] + matrix_transloc["u31"] * float(self.x) + matrix_transloc["u32"] * float(self.y) + \
                    matrix_transloc["u33"] * float(self.z)

        self.x = atomx
        self.y = atomy
        self.z = atomz



    def get_centre2atoms(self, atom):

        atomcentral = deepcopy(self)
        atomcentral.serial = self.serial
        atomcentral.name = "X"
        atomcentral.occupancy = 0.0
        atomcentral.Bfact = 0.0
        atomcentral.element = "X"
        atomcentral.x = (float(self.x) + float(atom.x)) /2
        atomcentral.y = (float(self.x) + float(atom.y)) / 2
        atomcentral.z = (float(self.x) + float(atom.z)) / 2

        return atomcentral

    def get_type_atom(self):

        try:
            return self.type
        except:

            if self.resName == "GLU" or self.resName == "ASP":
                if self.name == "OE1" or self.name == "OE2" or self.name == "OD1" or self.name == "OD2":
                    self.type = "Oox"
                    return self.type


            if self.resName == "TYR":
                if self.name == "OH":
                    self.type = "Oph"
                    return self.type

            if self.resName == "THR":
                if self.name == "OG1":
                    self.type = "Oh"
                    return self.type

            if self.resName == "SER":
                if self.name == "OG":
                    self.type = "Oh"
                    return self.type

            # Nitrogen histidine
            if self.resName == "HIS":
                if self.name == "NE2" or self.name  == "ND1":
                    self.type = "Nim"
                    return self.type

            # Nitrogen basic
            if self.resName == "LYS":
                if self.name == "NZ":
                    self.type = "NaI"
                    return self.type

            if self.resName == "ARG":
                if self.name == "NH1" or self.name == "NH2" or self.name == "NHE" or self.name == "NE":
                    self.type = "Ngu"
                    return self.type

                    #     if atom["resName"] == "ARG" :
                    #         if atom["name"] == "CZ" :
                    #             return "Cgu"

            # Nitrogen donor
            if self.resName == "ASN":
                if self.name == "ND2":
                    self.type = "Nam"
                    return self.type

            if self.resName == "GLN":
                if self.name == "NE2":
                    self.type = "Nam"
                    return self.type

            if self.resName in LRES:
                if self.name == "N":
                    self.type = "Nam"
                    return self.type

            # Carbon aromatic
            if self.resName == "PHE" or self.resName == "TYR":
                if self.name != "CA":
                    if self.name != "C":
                        self.type = "Car"
                        return self.type

            if self.resName == "TRP":
                if self.name != "CA":
                    if self.name != "CB":
                        self.type = "Car"
                        return self.type
            # O peptitique
            if self.resName in LRES:
                if self.name == "O":
                    self.type = "Oc"
                    return self.type

            # O acid amide
            if self.resName == " ASN" or self.resName == "GLN":
                if self.name == "OD1" or self.name == "OE1":
                    self.type = "Oc"
                    return self.type

            # water
            if self.resName == "HOH":
                if self.name == "O":
                    self.type = "Ow"
                    return self.type

            # C peptidic why ?
            #if self.resName in LRES:
            #    if self.name == "C":
            #        return "Oc"

            # sulfur
            if self.resName in LRES:
                if self.element == "S":
                    self.type = "Su"
                    return self.type

            self.type = "Xot"
            return self.type


    def builder_atom (self, name, serial, char, resName, chainID, resSeq, iCode, x, y, z, element, charge, occupancy, Bfact, recoder):

        self.name = name
        self.serial = serial
        self.char = char
        self.resName = resName
        self.chainID = chainID
        self.resSeq = resSeq
        self.iCode = iCode
        self.x = x
        self.y = y
        self.z = z
        # case of element not present in models
        if element != "":
            self.element = element
        else:
            self.element = name[0]
        #print self.element, element, "l487"
        self.charge = charge
        self.occupancy = occupancy
        self.Bfact = Bfact
        self.recoder = recoder


    def writeAtom(self, filout):

        filout.write("%-6s%5s %4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n" % (
            self.recoder, self.serial, self.name, self.char, self.resName, self.chainID, self.resSeq, self.iCode, self.x,
            self.y, self.z, self.occupancy, self.Bfact, self.element, self.charge))




class MD_PDB:

    def __init__(self, pPDB):
        self.pDM = pPDB

    def splitDM(self):
        prout = self.pDM.split(".")[0] + "/"
        self.lpdbfiles = []
        try: makedirs(prout)
        except: pass

        filin = open(self.pDM, "r")
        stread = filin.read()
        filin.close()

        lmodel = stread.split("MODEL")
        print len(lmodel)

        i = 2# remove first model
        imax = len(lmodel)
        #imax = 3 ############ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        while i <= imax:
            pfilout = prout + "MODEL-" + str(i-1) + ".pdb"
            self.lpdbfiles.append(pfilout)
            filout = open(pfilout, "w")
            lmodeli = lmodel[i-1].split("\n")
            filout.write("\n".join(lmodeli[1:]))
            filout.close()
            print i
            i += 1





class ProtIon:

    def __init__(self, PDBID, HOH):
        self.PDBID = PDBID
        self.HOH = HOH
        cPDB = PDB(PDBID)
        d_res = cPDB.get_byres()

        d_prot = {}
        for resID in d_res.keys():
            res = resID.split("_")[0]
            if res in LRES:
                d_prot[resID] = deepcopy(d_res[resID])
        self.atomLig = d_prot





class Ligand:
    def __init__(self, ligandID, PDBID):

        """

        :param ligandID: ligand PDB ID 3 characters
        :param PDBID: PDB ID, 4 characters
        """

        self.ligID = ligandID
        self.PDBID = PDBID

        cPDB = PDB(PDBID)
        d_res = cPDB.get_byres()

        d_lig = {}
        for resID in d_res.keys():
            if search(ligandID, resID):
                # control if it is not a amino acid modified -> distance less than 2.0A from atom
                if HETATOMinProtein(d_res, resID) == 0:
                    d_lig[resID] = deepcopy(d_res[resID])

        self.nbLig = len(d_lig.keys())
        self.atomLig = d_lig



def HETATOMinProtein(d_protein, resID):
    l_atomHET = d_protein[resID]
    for atomHET in l_atomHET:
        for res in d_protein.keys():
            # also case with metal to keep
            if res == resID or not res.split("_")[0] in LRES:
                continue
            else:
                for atomres in d_protein[res]:
                    dist = atomHET.euclidiendist(atomres)
                    if dist >= MAXDISTBYRES:
                        break
                    elif dist <= MAXCONNECT:
                        return 1
    return 0



def changeRecoder(latoms, recorder="HETATM"):

    for atom in latoms:
        atom.recoder = recorder



def convert_ListAtomtoList(latomin):
    lout = []

    for atomin in latomin:
        nameatom = str(atomin.resName) + str(atomin.resSeq)
        if not nameatom in lout:
            lout.append(nameatom)
    return lout


#dtest = PDB(PDB_input="/home/aborrel/MDPockets/dataWNK/pdb_WNK1_with_lig/WNK-out-MD/MODEL-1.pdb", hydrogen=1)
#dtest.change_bfactor(dp_bfactor="/home/aborrel/MDPockets/resultWNKLig/analysis/FPIBS/FPIScore.csv")
#dtest.writePDB(pfilout="/home/aborrel/MDPockets/resultWNKLig/analysis/FPIBS/protBfact.pdb")