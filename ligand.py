from pydpi.drug import constitution, topology, connectivity, kappa, bcut, basak, estate, moran, moreaubroto, geary, \
    charge, molproperty, moe, fingerprint
from pydpi import pydrug
from molvs import standardize_smiles, Standardizer


from rdkit.Chem.SaltRemover import SaltRemover
from rdkit import Chem

from copy import deepcopy
from os import path, getcwd, remove, system, listdir
from shutil import copy
from re import search
from numpy import mean, std

import toolbox
import pathFolder
import runExternalSoft
import descriptors3D



LSALT="[Co]"
LSMILESREMOVE=["[C-]#N", "[Al+3]", "[Gd+3]", "[Pt+2]", "[Au+3]", "[Bi+3]", "[Al]", "[Si+4]", "[Fe]", "[Zn]", "[Fe+2]",
               "[Ru+8]", "[Fe+]", "[Sr++]", "[Fe+3]", "[O--]", "[OH-]", "[Mn++]", "[La+3]", "[Lu+3]", "[SH-]", "[Pt+4]",
               "[Fe++]", "[W]", "[Cu+2]", "[Cr+3]", "[Tc+7]", "[Xe]", "[Tl+]", "[Zn+2]", "[F-]", "[C]", "[He]", "N#N",
               "O=O", "Cl[Ra]Cl", "[Mn+2]", "N#[N+][O-]", "II", "[Ga+3]", "[Mo+10]", "[Zn]", "[Fe]", "[Si+4]", "[Al]"]



LKAPA = ['kappa1', 'kappa2', 'kappa3', 'kappam1', 'kappam2', 'kappam3', 'phi']
LBUCUT =["bcutp16","bcutp15","bcutp14","bcutp13","bcutp12","bcutp11","bcutp10",
        "bcutp9","bcutp8","bcutp7","bcutp6","bcutp5","bcutp4","bcutp3",
        "bcutp2","bcutp1"]
LESTATE = ['Smax38', 'Smax39', 'Smax34', 'Smax35', 'Smax36', 'S43', 'Smax30', 'Smax31', 'Smax32', 'Smax33', 'S57',
           'S56', 'S55', 'S54', 'S53', 'S52', 'S51', 'S50', 'Smin49', 'S59', 'S58', 'Smin69', 'Smin68', 'Smin27',
           'Sfinger30', 'Sfinger31', 'Sfinger32', 'Sfinger33', 'Sfinger34', 'Sfinger35', 'Sfinger36', 'Sfinger37',
           'Sfinger38', 'Sfinger39', 'Smax2', 'Smax3', 'Smax4', 'Smax5', 'Smax6', 'Smax7', 'Smin77', 'Smax29', 'Smax37',
           'Smax23', 'Smax22', 'Smax21', 'Smax20', 'Smax27', 'Smax26', 'Smax25', 'Smax24', 'S44', 'S45', 'S46', 'S47',
           'S40', 'S41', 'S42', 'S17', 'Smin44', 'S48', 'S49', 'Smin8', 'Smin29', 'Smin28', 'Sfinger45', 'Sfinger44',
           'Sfinger47', 'Sfinger46', 'Sfinger41', 'Sfinger40', 'Sfinger43', 'Sfinger42', 'Smax47', 'Smin73', 'Smin70',
           'Smin71', 'Sfinger49', 'Sfinger48', 'Smin74', 'Smin75', 'Smin67', 'Smin6', 'Smin9', 'Smin7', 'Smin47',
           'Smax41', 'S79', 'S78', 'Smin19', 'Smax58', 'Smax59', 'S71', 'S70', 'S73', 'S72', 'S75', 'S74', 'S77',
           'S76', 'Smax73', 'Smin78', 'Sfinger56', 'Sfinger57', 'Sfinger54', 'Sfinger55', 'Sfinger52', 'Sfinger53',
           'Sfinger50', 'Sfinger51', 'Smin61', 'Smin60', 'Smin63', 'Smin62', 'Smin65', 'Smin64', 'Sfinger58',
           'Sfinger59', 'Smin48', 'Smin42', 'Smin76', 'Smin41', 'Smin72', 'Smax40', 'Smin40', 'Smax49', 'Smax48',
           'S68', 'S69', 'S66', 'S67', 'S64', 'S65', 'S62', 'S63', 'S60', 'S61', 'Smin54', 'Smax52', 'Sfinger69',
           'Sfinger68', 'Smin50', 'Smin51', 'Smin52', 'Smin53', 'Sfinger63', 'Sfinger62', 'Sfinger61', 'Sfinger60',
           'Sfinger67', 'S10', 'Sfinger65', 'Sfinger64', 'S13', 'S12', 'Sfinger76', 'Smin56', 'S9', 'S8', 'S3', 'S2',
           'S1', 'Smin55', 'S7', 'S6', 'S5', 'S4', 'Smax78', 'Smax45', 'Smax11', 'Sfinger72', 'Smin66', 'Smax44',
           'Smax70', 'Smax71', 'Smax72', 'S14', 'Smax74', 'Smax75', 'Smax76', 'Smax77', 'Smin43', 'Smax8', 'S19',
           'S18', 'Sfinger78', 'Sfinger79', 'Smin45', 'Smax9', 'Sfinger74', 'Sfinger75', 'S11', 'Sfinger77',
           'Sfinger70', 'Sfinger71', 'S15', 'Sfinger73', 'Smax43', 'Smin16', 'Smax42', 'Smax53', 'Smax66', 'Smax65',
           'Smax64', 'Smax63', 'Smax62', 'Smax61', 'Smax60', 'Smin26', 'Smax69', 'Smax68', 'Smax0', 'Smin57', 'Smax1',
           'Smin17', 'Smin36', 'Smin37', 'Smin34', 'Smin35', 'Smin32', 'Smin33', 'Smin30', 'Smin31', 'Smax67', 'Smin46',
           'Smax51', 'Smin38', 'Smin39', 'Smax12', 'Smax13', 'Smax10', 'S16', 'Smax16', 'Smax17', 'Smax14', 'Smax15',
           'Smin20', 'Smax18', 'Smax19', 'Sfinger66', 'Smax56', 'Smax28', 'Smax57', 'Smax54', 'Smin58', 'Smax55', 'S39',
           'S38', 'Smax46', 'S35', 'S34', 'S37', 'S36', 'S31', 'S30', 'S33', 'S32', 'Smin25', 'Smin24', 'Sfinger18',
           'Sfinger19', 'Smin21', 'Smax50', 'Smin23', 'Smin22', 'Sfinger12', 'Sfinger13', 'Sfinger10', 'Sfinger11',
           'Sfinger16', 'Sfinger17', 'Sfinger14', 'Sfinger15', 'Sfinger8', 'Sfinger9', 'Smin4', 'Smin5', 'Smin2',
           'Smin3', 'Smin0', 'Smin1', 'Sfinger1', 'Sfinger2', 'Sfinger3', 'Sfinger4', 'Sfinger5', 'Sfinger6',
           'Sfinger7', 'S22', 'S23', 'S20', 'S21', 'S26', 'S27', 'S24', 'S25', 'Smin59', 'S28', 'S29', 'Smin18',
           'Smin10', 'Smin11', 'Smin12', 'Smin13', 'Smin14', 'Smin15', 'Sfinger29', 'Sfinger28', 'Sfinger27',
           'Sfinger26', 'Sfinger25', 'Sfinger24', 'Sfinger23', 'Sfinger22', 'Sfinger21', 'Sfinger20']
LMOREAUBROTO = ['ATSe1', 'ATSe2', 'ATSe3', 'ATSe4', 'ATSe5', 'ATSe6', 'ATSe7', 'ATSe8', 'ATSp8', 'ATSp3', 'ATSv8',
                'ATSp1', 'ATSp7', 'ATSp6', 'ATSp5', 'ATSp4', 'ATSv1', 'ATSp2', 'ATSv3', 'ATSv2', 'ATSv5', 'ATSv4',
                'ATSv7', 'ATSv6', 'ATSm8', 'ATSm1', 'ATSm2', 'ATSm3', 'ATSm4', 'ATSm5', 'ATSm6', 'ATSm7']
LMORAN = ['MATSv8', 'MATSp4', 'MATSp8', 'MATSv1', 'MATSp6', 'MATSv3', 'MATSv2', 'MATSv5', 'MATSv4', 'MATSv7', 'MATSv6',
          'MATSm8', 'MATSp1', 'MATSm4', 'MATSm5', 'MATSm6', 'MATSm7', 'MATSm1', 'MATSm2', 'MATSm3', 'MATSe4', 'MATSe5',
          'MATSe6', 'MATSe7', 'MATSe1', 'MATSe2', 'MATSe3', 'MATSe8', 'MATSp3', 'MATSp7', 'MATSp5', 'MATSp2']
LGEARY = ['GATSp8', 'GATSv3', 'GATSv2', 'GATSv1', 'GATSp6', 'GATSv7', 'GATSv6', 'GATSv5', 'GATSv4', 'GATSe2', 'GATSe3',
          'GATSv8', 'GATSe6', 'GATSe7', 'GATSe4', 'GATSe5', 'GATSp5', 'GATSp4', 'GATSp7', 'GATSe1', 'GATSp1', 'GATSp3',
          'GATSp2', 'GATSe8', 'GATSm2', 'GATSm3', 'GATSm1', 'GATSm6', 'GATSm7', 'GATSm4', 'GATSm5', 'GATSm8']
LMOE = ['EstateVSA8', 'EstateVSA9', 'EstateVSA4', 'EstateVSA5', 'EstateVSA6', 'EstateVSA7', 'EstateVSA0', 'EstateVSA1',
        'EstateVSA2', 'EstateVSA3', 'PEOEVSA13', 'PEOEVSA12', 'PEOEVSA11', 'PEOEVSA10', 'MTPSA', 'VSAEstate0',
        'VSAEstate1', 'VSAEstate2', 'VSAEstate3', 'VSAEstate4', 'VSAEstate5', 'VSAEstate6', 'VSAEstate7', 'VSAEstate8',
        'LabuteASA', 'PEOEVSA3', 'PEOEVSA2', 'PEOEVSA1', 'PEOEVSA0', 'PEOEVSA7', 'PEOEVSA6', 'PEOEVSA5', 'PEOEVSA4',
        'MRVSA5', 'MRVSA4', 'PEOEVSA9', 'PEOEVSA8', 'MRVSA1', 'MRVSA0', 'MRVSA3', 'MRVSA2', 'MRVSA9', 'slogPVSA10',
        'slogPVSA11', 'MRVSA8', 'MRVSA7', 'MRVSA6', 'EstateVSA10', 'slogPVSA2', 'slogPVSA3', 'slogPVSA0', 'slogPVSA1',
        'slogPVSA6', 'slogPVSA7', 'slogPVSA4', 'slogPVSA5', 'slogPVSA8', 'slogPVSA9', 'VSAEstate9', 'VSAEstate10']



class Descriptors:
    def __init__(self, dcompound, logfile, writecheck=1, kSMILES="CANONICAL_SMILES", kID="CMPD_CHEMBLID"):
        self.compound = dcompound
        self.kID = kID
        loader = pydrug.PyDrug()

        # if SMILES, load using SMILES code
        if not kSMILES in dcompound.keys():
            try:
                smile = runExternalSoft.babelConvertSDFtoSMILE(dcompound["sdf"])
                self.compound[kSMILES] = smile
            except:
                print "ERROR INPUT SDF - l33"
                self.log = "ERROR"
                try:logfile.write(
                    self.compound[kID] + "\t---\tERROR-SDF ORIGINAL INPUT\n")
                except:pass
                return


        #Standardize smile code
        try: smilestandadized = standardize_smiles(self.compound[kSMILES])
        except:
            logfile.write(self.compound[kID] + "\t" + str(self.compound[kSMILES]) + "\tERROR-SMILES INPUT"
                                                                                               "\n")
            self.log = "ERROR"
            return

        #Standardize using molvs (http://molvs.readthedocs.io/en/latest/api.html#molvs-fragment)
        s = Standardizer()
        mol = Chem.MolFromSmiles(smilestandadized)
        molstandardized = s.standardize(mol)
        smilestandadized = Chem.MolToSmiles(molstandardized)

        # remove salt
        # 1.default
        remover = SaltRemover()
        mol = Chem.MolFromSmiles(smilestandadized)
        molcleandefault = remover(mol)
        # 2. Personal remover
        homeremover = SaltRemover(defnData=LSALT)
        molclean = homeremover(molcleandefault)
        smilesclean = Chem.MolToSmiles(molclean)
        # 3. SMILES remove other manual salts + fragments -> for fragment take one if exactly same compound
        lelem = smilesclean.split(".")
        if len(lelem) > 1:
            # reduce double, case of several salts are included - 255
            lelem = list(set(lelem))
            for smilesdel in LSMILESREMOVE:
                if smilesdel in lelem:
                    lelem.remove(smilesdel)
            try:lelem.remove("")# case of bad smile
            except:pass
            if len(lelem) == 1:
                smilesclean = str(lelem[0])
            else:
                # 4. Fragments
                #Case of fragment -> stock in log file, check after to control
                logfile.write(self.compound[kID] + "\t" + str(self.compound[kSMILES]) + "\tFRAGMENT IN INPUT"
                                                                                                   "\n")
                print ".".join(lelem), " - FRAGMENTS - l66"
                self.log = "ERROR"
                return
        else:
            pass


        print self.compound[kSMILES], "SMILES IN - l25 liganddescriptors"
        print smilesclean, "SMILES without salt and standardized"

        # case where only salt are included
        if smilesclean == "":
            logfile.write(self.compound[kID] + "\t" + str(self.compound[kSMILES]) + "\tEMPTY SMILES AFTER "
                                                                                               "STANDARDIZATION\n")
            print "EMPTY SMILES AFTER STANDARDIZATION - l84"
            self.log = "ERROR"
            return

        self.compound[kSMILES] = smilesclean
        self.log = "OK"

        if writecheck == 1:
            # SMILES code
            pfileSMILES = pathFolder.PR_COMPOUNDS + str(dcompound[kID]) + ".smi"
            fileSMILES = open(pfileSMILES, "w")
            fileSMILES.write(self.compound[kSMILES])
            fileSMILES.close()

            # SDF input
            if "sdf" in self.compound.keys():
                pfileSDF = pathFolder.PR_COMPOUNDS + str(dcompound[kID]) + ".sdf"
                fileSDF = open(pfileSDF, "w")
                fileSDF.write(self.compound["sdf"])
                fileSDF.close()

        # read mol
        self.mol = loader.ReadMolFromSmile(self.compound[kSMILES])


    def get_descriptor1D2D(self):

        try:
            self.consti = constitution.GetConstitutional(self.mol)
        except:
            self.consti = {}
        self.compo = {}
        try:
            self.compo["nheavy"] = self.mol.GetNumHeavyAtoms()
        except:
            self.compo = {}
        try:
            self.molprop = molproperty.GetMolecularProperty(self.mol)
        except:
            self.molprop = {}
        try:
            self.topo = topology.GetTopology(self.mol)
        except:
            self.topo = {}
        try:
            self.connect = connectivity.GetConnectivity(self.mol)
        except:
            self.connect = {}
        try:
            self.kap = kappa.GetKappa(self.mol)
        except:
            self.kap = {}
        try:
            self.burden = bcut.GetBurden(self.mol)
        except:
            self.burden = {}
        try:
            self.basakD = basak.Getbasak(self.mol)
        except:
            self.basakD = {}
        try:
            self.est = estate.GetEstate(self.mol)
        except:
            self.est = {}
        try:
            self.moreauBurto = moreaubroto.GetMoreauBrotoAuto(self.mol)
        except:
            self.moreauBurto = {}
        try:
            self.autcormoran = moran.GetMoranAuto(self.mol)
        except:
            self.autcormoran = {}
        try:
            self.gearycor = geary.GetGearyAuto(self.mol)
        except:
            self.gearycor = {}
        try:
            self.charges = charge.GetCharge(self.mol)
        except:
            self.charges = {}
        try:
            self.MOE = moe.GetMOE(self.mol)
        except:
            self.MOE = {}

        # list 1D2D -> modified in main library !!!!
        self.l1D2D = constitution._constitutional.keys()
        self.l1D2D = self.l1D2D + ["nheavy"]
        self.l1D2D = self.l1D2D + molproperty.MolecularProperty.keys()
        self.l1D2D = self.l1D2D +topology._Topology.keys()
        self.l1D2D = self.l1D2D + connectivity._connectivity.keys()
        self.l1D2D = self.l1D2D + LKAPA
        self.l1D2D = self.l1D2D + LBUCUT
        self.l1D2D = self.l1D2D + basak._basak.keys()
        self.l1D2D = self.l1D2D + LESTATE
        self.l1D2D = self.l1D2D + LMOREAUBROTO
        self.l1D2D = self.l1D2D + LMORAN
        self.l1D2D = self.l1D2D + LGEARY
        self.l1D2D = self.l1D2D + charge._Charge.keys()
        self.l1D2D = self.l1D2D + LMOE

        # combine all 1D-2D
        self.all1D2D = dict()
        self.all1D2D.update(self.consti)
        self.all1D2D.update(self.compo)
        self.all1D2D.update(self.molprop)
        self.all1D2D.update(deepcopy(self.topo))
        self.all1D2D.update(deepcopy(self.connect))
        self.all1D2D.update(deepcopy(self.kap))
        self.all1D2D.update(deepcopy(self.burden))
        self.all1D2D.update(deepcopy(self.basakD))
        self.all1D2D.update(deepcopy(self.est))
        self.all1D2D.update(deepcopy(self.moreauBurto))
        self.all1D2D.update(deepcopy(self.autcormoran))
        self.all1D2D.update(deepcopy(self.gearycor))
        self.all1D2D.update(deepcopy(self.charges))
        self.all1D2D.update(deepcopy(self.MOE))


    def get_descriptor3D(self, pr3D):

        #if pr3D != "":
        #    self.all3D = {}
        #else:
            # extract docking poses
            CHEMBLID = self.compound[self.kID]
            psdf3D = pr3D + CHEMBLID + ".1.sdf"

            self.all3D = descriptors3D.get3Ddesc(psdf3D)
            self.l3D = descriptors3D.l3D




    def writeTablesDesc(self, presult, typeDesc, kSMILES="CANONICAL_SMILES", kID="CMPD_CHEMBLID", unique=0):


        # case we would like one unique file, if extention is .csv
        if "all1D2D" in self.__dict__ and typeDesc == "1D2D":
            ldesc = self.l1D2D
        elif "all3D" in self.__dict__ and typeDesc == "3D":
            ldesc = self.l3D
        else:
            ldesc = self.l1D2D + self.l3D

        if not path.exists(presult):
            filout = open(presult, "w")
            # header
            filout.write("ID\tSMILES\t")

            filout.write("\t".join(ldesc) + "\n")

        else:
            filout = open(presult, "a")

        filout.write(self.compound[kID] + "\t" + self.compound[kSMILES])

        for desc in ldesc:
            try:
                filout.write("\t" + str(self.all1D2D[desc]))
            except:
                filout.write("\t" + str(self.all3D[desc]))

        filout.write('\n')
        filout.close()


