from os import path,listdir
from re import search

import runExternalSoft
import toolbox
import pathFolder




class MDQSAR():

    def __init__(self, prMDdescc, prout):
        self.prDesc = prMDdescc
        self.prout = prout


    def builtDataset(self):

        # add control exist or not
        #

        ddesc = {}

        lprDMdesc = listdir(self.prDesc)
        for prMDdesc in lprDMdesc:
            MDID = prMDdesc
            lfilesDesc = listdir(self.prDesc + prMDdesc + "/")
            for fileDesc in lfilesDesc:
                if search("^desc", fileDesc):
                    typeDesc = fileDesc[4:]
                    if not typeDesc in ddesc.keys():
                        ddesc[typeDesc] = {}
                    ddescMIID = toolbox.loadTable(self.prDesc + prMDdesc + "/" + fileDesc)
                    for desc in ddescMIID[MDID].keys():
                        if not desc in ddesc[typeDesc].keys():
                            ddesc[typeDesc][desc] = {}
                        ddesc[typeDesc][desc][MDID] = ddescMIID[MDID][desc]

        self.dataset = ddesc



    def writeDataset(self,ltypeDesc, prout):

        if not "dataset" in self.__dict__:
            self.builtDataset()

        pfilout = prout + "desc" + "-".join(ltypeDesc)
        filout = open(pfilout, "w")

        #ldesc
        ldescWrite = []
        lMDID = []
        for typeDesc in ltypeDesc:
            ldescWrite = ldescWrite + self.dataset[typeDesc].keys()
            for desc in self.dataset[typeDesc].keys():
                for MDID in self.dataset[typeDesc][desc].keys():
                    if not MDID in lMDID:
                        lMDID.append(MDID)


        filout.write("MDID\t" + "\t".join(ldescWrite) + "\n")

        for MDID in lMDID:
            filout.write(str(MDID))
            for desc in ldescWrite:
                flag = 0
                for typedesc in ltypeDesc:
                    try:
                        filout.write("\t" + str(self.dataset[typedesc][desc][MDID]))
                        flag = 1
                        break
                    except:
                        pass

                if flag == 0 :
                    filout.write("\tNA")
            filout.write("\n")
        filout.close()




    def datasetAnalysis(self):

        return


    def runQSARModel(self, ltypeDesc, paffinity, corcoef, maxQuantile, valSplit):

        prout = self.prout + "-".join(ltypeDesc) + "/"
        pathFolder.createFolder(prout)

        pdesc = self.writeDataset(ltypeDesc, prout)

        dfileTrainTest = runExternalSoft.prepareDataset(pdesc, paffinity_currated, prQSAR, corcoef=corcoef,
                                                        maxQuantile=maxQuantile, valSplit=valSplit)

        #runExternalSoft.



        return


