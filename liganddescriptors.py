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

import toolbox
import pathFolder
import runExternalSoft

LSALT="[Co]"
LSMILESREMOVE=["[C-]#N", "[Al+3]", "[Gd+3]", "[Pt+2]", "[Au+3]", "[Bi+3]", "[Al]", "[Si+4]", "[Fe]", "[Zn]", "[Fe+2]",
               "[Ru+8]", "[Fe+]", "[Sr++]", "[Fe+3]", "[O--]", "[OH-]", "[Mn++]", "[La+3]", "[Lu+3]", "[SH-]", "[Pt+4]",
               "[Fe++]", "[W]", "[Cu+2]", "[Cr+3]", "[Tc+7]", "[Xe]", "[Tl+]", "[Zn+2]", "[F-]", "[C]", "[He]", "N#N",
               "O=O", "Cl[Ra]Cl", "[Mn+2]", "N#[N+][O-]", "II", "[Ga+3]", "[Mo+10]", "[Zn]", "[Fe]", "[Si+4]", "[Al]"]

class Descriptors:
    def __init__(self, dcompound, logfile, writecheck=1, kSMILES="CANONICAL_SMILES", kID="CMPD_CHEMBLID"):
        self.compound = dcompound
        loader = pydrug.PyDrug()

        # if SMILES, load using SMILES code
        if not kSMILES in dcompound.keys():
            try:
                smile = toolbox.babelConvertSDFtoSMILE(dcompound["sdf"])
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

    def get_descriptorOD1D(self):
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

        # combine all 1D
        self.all1D = {}
        self.all1D.update(self.consti)
        self.all1D.update(self.compo)
        self.all1D.update(self.molprop)

        # listdesc
        self.l1D = constitution._constitutional.keys()
        self.l1D = self.l1D + ["nheavy"]
        self.l1D = self.l1D + molproperty.MolecularProperty.keys()

    def get_descriptor2D(self):
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

        # list 2D -> modified in main library !!!!
        self.l2D = topology._Topology.keys()
        self.l2D = self.l2D + connectivity._connectivity.keys()
        self.l2D = self.l2D + kappa._kapa.keys()
        self.l2D = self.l2D + bcut._bcut
        self.l2D = self.l2D + basak._basak.keys()
        self.l2D = self.l2D + estate._estate.keys()
        self.l2D = self.l2D + moreaubroto._moreaubroto.keys()
        self.l2D = self.l2D + moran._moran.keys()
        self.l2D = self.l2D + geary._geary.keys()
        self.l2D = self.l2D + charge._Charge.keys()
        self.l2D = self.l2D + moe._moe.keys()

        # combine all 2D
        self.all2D = dict()
        self.all2D.update(deepcopy(self.topo))
        self.all2D.update(deepcopy(self.connect))
        self.all2D.update(deepcopy(self.kap))
        self.all2D.update(deepcopy(self.burden))
        self.all2D.update(deepcopy(self.basakD))
        self.all2D.update(deepcopy(self.est))
        self.all2D.update(deepcopy(self.moreauBurto))
        self.all2D.update(deepcopy(self.autcormoran))
        self.all2D.update(deepcopy(self.gearycor))
        self.all2D.update(deepcopy(self.charges))
        self.all2D.update(deepcopy(self.MOE))

    def get_fingerprints(self):
        # fingerprint
        self.fingerAtomPairs = fingerprint.CalculateAtomPairsFingerprint(self.mol)
        self.fingerDaylight = fingerprint.CalculateDaylightFingerprint(self.mol)
        self.fingerEstate = fingerprint.CalculateEstateFingerprint(self.mol)
        self.fingerFP4 = fingerprint.CalculateFP4Fingerprint(self.mol)
        self.fingerMACCS = fingerprint.CalculateMACCSFingerprint(self.mol)
        self.fingerMorgan = fingerprint.CalculateMorganFingerprint(self.mol)
        self.fingerTorsion = fingerprint.CalculateTopologicalTorsionFingerprint(self.mol)

    def get_descriptor3D(self, log):
        """
        Compute descriptors 3D from SMILES code and generate the 3D using ligprep
        :return: dictionary of descriptors in all3D
        """

        # clean temp folder - used to compute 3D descriptors
        prtemp = pathFolder.cleanFolder()
        psdf3Dout = pathFolder.PR_COMPOUNDS + str(self.compound["DATABASE_ID"]) + "_3D.sdf"

        # temp SMILES
        pfilesmile = prtemp + "tem.smi"
        filesmile = open(pfilesmile, "w")
        filesmile.write(self.compound["SMILES"])
        filesmile.close()

        # run ligprep
        if not path.exists(psdf3Dout):
            psdf3D = runExternalSoft.runLigprep(psmilin=pfilesmile)

            # case error in ligprep
            if not path.exists(psdf3D) or path.getsize(psdf3D) == 0:
                self.all3D = toolbox.parsePadelOut()
                self.l3D = self.all3D.keys()
                log.write(self.compound["DATABASE_ID"] + "\t" + self.compound["SMILES"] + "\t" + psdf3D)
                pdesc = ""
            else:
                psdf3Dout = toolbox.selectMinimalEnergyLigPrep(psdfin=psdf3D,
                                                               psdfout=psdf3Dout)
                # take only best energy
                remove(psdf3D)
                remove(pfilesmile)
                copy(psdf3Dout, psdf3D)
                pdesc = runExternalSoft.runPadel(prtemp)

        # run 3D descriptor using Padel
        self.all3D = toolbox.parsePadelOut(pdesc)
        self.l3D = self.all3D.keys()


    def writeTablesDesc(self, prresult, kSMILES="CANONICAL_SMILES", kID="CMPD_CHEMBLID", unique=0):


        # case we would like one unique file, if extention is .csv
        if prresult[-3:] == "csv":
            if not path.exists(prresult):
                filout = open(prresult, "w")
                # header
                filout.write("ID\tSMILES\t")
                lheader = []
                if "all1D" in self.__dict__:
                    lheader = lheader + self.l1D
                if "all2D" in self.__dict__:
                    lheader = lheader + self.l2D
                if "all3D" in self.__dict__:
                    lheader = lheader + self.l3D
                filout.write("\t".join(lheader) + "\n")
            else:
                filout = open(prresult, "a")

            filout.write(self.compound[kID] + "\t" + self.compound[kSMILES])

            if "all1D" in self.__dict__:
                for desc1D in self.l1D:
                    try:
                        filout.write("\t" + str(self.all1D[desc1D]))
                    except:
                        filout.write("\tNA")

            if "all2D" in self.__dict__:
                for desc2D in self.l2D:
                    try:
                        filout.write("\t" + str(self.all2D[desc2D]))
                    except:
                        filout.write("\tNA")

            if "all3D" in self.__dict__:
                for desc3D in self.l3D:
                    try:
                        filout.write("\t" + str(self.all2D[desc3D]))
                    except:
                        filout.write("\tNA")

            filout.write('\n')
            filout.close()

        else:
            # case 1D
            if "all1D" in self.__dict__:
                if not path.exists(prresult + "1D.csv"):
                    self.fil1D = open(prresult + "1D.csv", "w")
                    # header
                    self.fil1D.write("ID\tSMILES\t")
                    self.fil1D.write("\t".join(self.l1D) + "\n")
                else:
                    self.fil1D = open(prresult + "1D.csv", "a")
                self.fil1D.write(self.compound[kID] + "\t" +self.compound[kSMILES])

                for desc1D in self.l1D:
                    try:
                        self.fil1D.write("\t" + str(self.all1D[desc1D]))
                    except:
                        self.fil1D.write("\tNA")
                self.fil1D.write("\n")
                self.fil1D.close()

            # case 2D
            if "all2D" in self.__dict__:
                if not path.exists(prresult + "2D.csv"):
                    self.fil2D = open(prresult + "2D.csv", "w")
                    # header
                    self.fil2D.write("ID\tSMILES\t")
                    self.fil2D.write("\t".join(self.l2D) + "\n")
                else:
                    self.fil2D = open(prresult + "2D.csv", "a")

                self.fil2D.write(self.compound[kID] + "\t" + self.compound[kSMILES])
                for desc2D in self.l2D:
                    try:
                        self.fil2D.write("\t" + str(self.all2D[desc2D]))
                    except:
                        self.fil2D.write("\tNA")
                self.fil2D.write("\n")
                self.fil2D.close()

            # case 3D - not work
            if "all3D" in self.__dict__:
                if not path.exists(prresult + "3D.csv", ):
                    self.fil3D = open(prresult + "3D.csv", "w")
                    # header
                    self.fil3D.write("ID\tSMILES\t")
                    self.fil3D.write("\t".join(self.l3D) + "\n")
                else:
                    self.fil3D = open(prresult + "3D.csv", "a")

                self.fil3D.write(self.compound[kID] + "\t" + self.compound[kSMILES])
                for desc3D in self.l3D:
                    try:
                        self.fil3D.write("\t" + str(self.all3D[desc3D]))
                    except:
                        self.fil3D.write("\tNA")
                self.fil3D.write("\n")
                self.fil3D.close()
