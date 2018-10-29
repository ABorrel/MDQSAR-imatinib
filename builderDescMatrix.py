from os import path, listdir
from re import search

import runExternalSoft
import toolbox
import pathFolder





class Builder:
    def __init__(self, prMDdescc, pdesc2D, pdesc3D, paffinity, corcoef, maxquantile, valsplit, typeAff):
        self.prMDDesc = prMDdescc
        self.paff = paffinity
        self.pdesc2D = pdesc2D
        self.pdesc3D = pdesc3D
        self.typeAff = typeAff

        # parameters
        self.corcoef = corcoef
        self.maxQuantile = maxquantile
        self.valsplit = valsplit


    def buildDataset(self, ltypeDesc, typeAff, prout):

        # descriptor to select
        self.ltypedesc = ltypeDesc
        self.prout = prout
        self.typeAff = typeAff

        ddesc = {}
        for typeDesc in self.ltypedesc:
            if typeDesc == "Lig2D":
                ddesLig = toolbox.loadTable(self.pdesc2D)
                ddesc.update(ddesLig)
                #print ddesc.keys()[0]
                #print ddesc[ddesc.keys()[0]].keys()[1:10]
            elif typeDesc == "Lig3D":
                ddesLig = toolbox.loadTable(self.pdesc3D)
                ddesc.update(ddesLig)
            else:
                lprDMCompound = listdir(self.prMDDesc)
                for prDMCompound in lprDMCompound:
                    lfilesDesc = listdir(self.prMDDesc + prDMCompound + "/")
                    for fileDesc in lfilesDesc:
                        if search("^desc", fileDesc) and search(typeDesc, fileDesc):
                            ddesctemp = toolbox.loadTable(self.prMDDesc + prDMCompound + "/" + fileDesc)
                            ddesc.update(ddesctemp)

        #print ddesc.keys()
        #print len(ddesc.keys())
        #print ddesc[ddesc.keys()[0]].keys()
        #print len(ddesc[ddesc.keys()[0]].keys())
        self.dataset = ddesc




    def writeDataset(self):

        if not "dataset" in self.__dict__:
            self.buildDataset()

        pfilout = self.prout + "descGlobal"
        # shortcut regeneration data
        if path.exists(pfilout) and path.getsize(pfilout) > 50:
            self.pdescglobal = pfilout
            return self.pdescglobal

        filout = open(pfilout, "w")

        #extract list descriptors
        ldesc = []
        for ID in self.dataset.keys():
            if len(self.dataset[ID].keys()) > len(ldesc):
                ldesc = self.dataset[ID].keys()

        # remove SMILES column
        if "SMILES" in ldesc:
            del ldesc[ldesc.index("SMILES")]
            #ldesc = ["SMILES"] + ldesc

        # header
        filout.write("ID\t" + "\t".join(ldesc) + "\n")

        for ID in self.dataset.keys():
            filout.write(str(ID))
            for desc in ldesc:
                try: filout.write("\t" + str(self.dataset[ID][desc]))
                except: filout.write("\tNA")
            filout.write("\n")
        filout.close()

        self.pdescglobal = pfilout
        return pfilout



    def prepMatrixDesc(self, pdesc, prout):

        pdesc = runExternalSoft.prepareMatrixDesc(pdesc, corcoef=self.corcoef, maxQuantile=self.maxQuantile, prout=prout)
        return pdesc


    def prepTrainTestSet(self):

        if not "dtrain" in self.__dict__:
            dfileTrainTest = runExternalSoft.prepareDataset(self.pdescglobal, self.paff, self.prout, corcoef=self.corcoef,
                                                    maxQuantile=self.maxQuantile, valSplit=self.valsplit,
                                                    typeAff=self.typeAff)

            print dfileTrainTest
            self.dtrain = dfileTrainTest["train"]
            self.dtest = dfileTrainTest["test"]

        else:
            dfileTrainTest = {}
            dtrain = runExternalSoft.createSetFromTable(self.pdescglobal, self.dtrain, self.prout + "trainSet_global.csv", self.prout)
            dfileTrain = runExternalSoft.prepareDataset(dtrain, self.paff, self.prout,
                                                            corcoef=self.corcoef,
                                                            maxQuantile=self.maxQuantile, valSplit=1,
                                                            typeAff=self.typeAff)
            dfileTrainTest["train"] = dfileTrain["train"]

            dtest = runExternalSoft.createSetFromTable(self.pdescglobal, self.dtest, self.prout + "testSet_global.csv", self.prout)
            dfileTest = runExternalSoft.prepareDataset(dtest, self.paff, self.prout,
                                                            corcoef=self.corcoef,
                                                            maxQuantile=self.maxQuantile, valSplit=1,
                                                            typeAff=self.typeAff)
            dfileTrainTest["test"] = dfileTest["train"]

        return dfileTrainTest


    def datasetAnalysis(self):
        """visualisation"""

        # create folder for visu
        prout = pathFolder.createFolder(self.prout + "VisuData/")
        runExternalSoft.QSARsVisuData(self.pdescglobal, self.paff, prout, self.corcoef, self.maxQuantile, "0")

        return
