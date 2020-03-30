import runExternalSoft
import pathFolder





class QSARModeling:
    def __init__(self, ltypeDesc, typeAff, prout):


        self.ltypedesc = ltypeDesc

        # define folder with QSAR
        prout = prout + "-".join(self.ltypedesc) + "_" + str(typeAff) + "/"
        pathFolder.createFolder(prout)
        self.prout = prout
        self.typeAff = typeAff


    def prep(self, builder):

        builder.buildDataset(self.ltypedesc, self.typeAff, self.prout)
        builder.writeDataset()
        dfileTrainTest = builder.prepTrainTestSet()
        self.dtrainTest = dfileTrainTest



    def runQSARModel(self):

        runExternalSoft.QSARsReg(self.dtrainTest["train"], self.dtrainTest["test"], "0", self.prout)


