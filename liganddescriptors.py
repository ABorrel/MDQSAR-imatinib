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
        self.l1D2D = self.l1D2D + kappa._kapa.keys()
        self.l1D2D = self.l1D2D + bcut._bcut
        self.l1D2D = self.l1D2D + basak._basak.keys()
        self.l1D2D = self.l1D2D + estate._estate.keys()
        self.l1D2D = self.l1D2D + moreaubroto._moreaubroto.keys()
        self.l1D2D = self.l1D2D + moran._moran.keys()
        self.l1D2D = self.l1D2D + geary._geary.keys()
        self.l1D2D = self.l1D2D + charge._Charge.keys()
        self.l1D2D = self.l1D2D + moe._moe.keys()

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

    def get_fingerprints(self):
        # fingerprint based on rdkit
        self.fingerAtomPairs = fingerprint.CalculateAtomPairsFingerprint(self.mol)
        self.fingerDaylight = fingerprint.CalculateDaylightFingerprint(self.mol)
        self.fingerEstate = fingerprint.CalculateEstateFingerprint(self.mol)
        self.fingerFP4 = fingerprint.CalculateFP4Fingerprint(self.mol)
        self.fingerMACCS = fingerprint.CalculateMACCSFingerprint(self.mol)
        self.fingerMorgan = fingerprint.CalculateMorganFingerprint(self.mol)
        self.fingerTorsion = fingerprint.CalculateTopologicalTorsionFingerprint(self.mol)

    def get_descriptor3D(self, lcpd, pdesc, pr3D=""):
        """
        Compute descriptors 3D from SMILES code and generate the 3D using ligprep
        :return: dictionary of descriptors in all3D
        """

        # from pose directory
        if pr3D != "":
            # extract docking poses
            lsdf = []
            prtemp = pathFolder.analyses(psub="desc/temp3D")
            for cpd in lcpd:
                psdf3D = pr3D + str(cpd["CMPD_CHEMBLID"]) + ".1.sdf"
                lsdf.append(str(cpd["CMPD_CHEMBLID"]) + ".1.sdf")
                print psdf3D, "l.239 - ligand descriptors"
                if path.exists(psdf3D):
                    psdftemp = prtemp + psdf3D.split("/")[-1]

                    #format sdf for MDLV3000Reader
                    cmdbabel = "babel " + psdf3D + " " + prtemp + str(cpd["CMPD_CHEMBLID"]) + ".1.sdf 2>/dev/null"
                    system(cmdbabel)

            # prepare par file
            pdesc = runExternalSoft.runKrakenX(prtemp, lsdf, pdesc)


        #else:
            # need to prepare 3D
            # clean temp folder - used to compute 3D descriptors
            #prtemp = pathFolder.cleanFolder()
            #psdf3Dout = pathFolder.PR_COMPOUNDS + str(self.compound["DATABASE_ID"]) + "_3D.sdf"

            # temp SMILES
            #pfilesmile = prtemp + "tem.smi"
            #filesmile = open(pfilesmile, "w")
            #filesmile.write(self.compound["SMILES"])
            #filesmile.close()

            # run ligprep - to develop
            #if not path.exists(psdf3Dout):
                #psdf3D = runExternalSoft.runLigprep(psmilin=pfilesmile)

                # case error in ligprep
                #if not path.exists(psdf3D) or path.getsize(psdf3D) == 0:
                #    pdesc = ""
                #else:
                #    psdf3Dout = toolbox.selectMinimalEnergyLigPrep(psdfin=psdf3D,
                #                                                       psdfout=psdf3Dout)
                    # take only best energy
                #    remove(psdf3D)
                #    remove(pfilesmile)
                #    copy(psdf3Dout, psdf3D)
                    #pdesc = runExternalSoft.runKrakenX(prtemp)

        # run 3D descriptor using Padel
        #self.all3D = toolbox.parseKrakenX(pdesc)
        #self.l3D = self.all3D.keys()


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
                    lheader = lheader + self.l1D2D
                if "all3D" in self.__dict__:
                    lheader = lheader + self.l3D
                filout.write("\t".join(lheader) + "\n")
            else:
                filout = open(prresult, "a")

            filout.write(self.compound[kID] + "\t" + self.compound[kSMILES])

            if "all1D" in self.__dict__:
                for desc1D in self.l1D:
                    try:
                        filout.write("\t" + str(self.all1D2D[desc1D]))
                    except:
                        filout.write("\tNA")

            if "all2D" in self.__dict__:
                for desc2D in self.l1D2D:
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
                        self.fil1D.write("\t" + str(self.all1D2D[desc1D]))
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
                    self.fil2D.write("\t".join(self.l1D2D) + "\n")
                else:
                    self.fil2D = open(prresult + "2D.csv", "a")

                self.fil2D.write(self.compound[kID] + "\t" + self.compound[kSMILES])
                for desc2D in self.l1D2D:
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



def MolecularDesc(lcpd, pfiloutdesc, prdocking = "", D1D = 1, D3D = 1, plog = "log.txt"):

    logfile = open(plog, "w")

    pdesc2D = pfiloutdesc + "2D.csv"
    if not path.exists(pdesc2D) and not path.getsize(pdesc2D) > 50:
        for compound in lcpd:
            dcompound = Descriptors(compound, logfile)

            if dcompound.log == "ERROR":
                continue

            if D1D == 1:
                dcompound.get_descriptor1D2D()

            dcompound.writeTablesDesc(pdesc2D)

    if D3D == 1:
        pdesc3D = pfiloutdesc + "3D.csv"
        if not path.exists(pdesc3D) and not path.getsize(pdesc3D) > 50:
            dcompound.get_descriptor3D(lcpd, prdocking, pdesc3D)
            # format table


    # combine Table
    if D3D == 1:
        d3D = toolbox.parseKrakenX(pdesc3D)
        print d3D.keys(), "CCCCC"
        ddesc = toolbox.loadTable(pdesc2D, d3D)
    else:
        ddesc = toolbox.loadTable(pdesc2D)

    print ddesc.keys()
    toolbox.writeTableDesc(ddesc, pfiloutdesc + "all.desc")

