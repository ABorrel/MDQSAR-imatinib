from os import path, listdir
from re import search

import runExternalSoft
import toolbox
import pathFolder





class QSARModeling:
    def __init__(self, prMDdescc, pdesc, paffinity, ltypeDesc, prout):
        self.prMDDesc = prMDdescc
        self.ltypedesc = ltypeDesc
        self.paff = paffinity
        self.pdesc = pdesc

        # define folder with QSAR
        prout = prout + "-".join(self.ltypedesc) + "/"
        pathFolder.createFolder(prout)
        self.prout = prout


    def builtDataset(self):


        ddesc = {}
        for typeDesc in self.ltypedesc:
            if typeDesc == "Lig2D":
                ddesLig = toolbox.loadTable(self.pdesc)
                ddesc.update(ddesLig)
                #print ddesc.keys()[0]
                #print ddesc[ddesc.keys()[0]].keys()[1:10]
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
            self.builtDataset()

        pfilout = self.prout + "descGlobal"
        # shortcut regeneration data
        #if path.exists(pfilout) and path.getsize(pfilout) > 50:
        #    self.pdescglobal = pfilout
        #    return self.pdescglobal

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



    def datasetAnalysis(self):

        "visualisation"

        return



    def runQSARModel(self, corcoef, maxquantile, valsplit):

        self.corcoef = corcoef
        self.maxQuantile = maxquantile
        self.valsplit = valsplit

        dfileTrainTest = runExternalSoft.prepareDataset(self.pdescglobal, self.paff, self.prout, corcoef=self.corcoef,
                                                        maxQuantile=self.maxQuantile, valSplit=self.valsplit)

        print dfileTrainTest
        runExternalSoft.QSARsReg(dfileTrainTest["train"], dfileTrainTest["test"], "0", self.prout)


