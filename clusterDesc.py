import pathFolder
import runExternalSoft
import toolbox

from os import path
from shutil import copyfile

class clusterDesc:


    def __init__(self, ltypeDesc, typeAff, cutoff, prout):

        self.ltypedesc = ltypeDesc
        self.cutoff = cutoff
        self.typeAff = typeAff

        # define folder with QSAR
        prout = prout + "-".join(self.ltypedesc) + "_" + str(typeAff) + "/"
        pathFolder.createFolder(prout)
        self.prout = prout




    def prep(self, builder):

        builder.buildDataset(self.ltypedesc, self.typeAff, self.prout)
        pdesc = builder.writeDataset()
        pdescClean = builder.prepMatrixDesc(pdesc, self.prout)

        self.paff = builder.paff
        self.pdesc = pdescClean


    def clusterize(self):

        runExternalSoft.clusterize(self.pdesc, self.paff, self.typeAff, self.cutoff, self.prout)


    def clusterizeTopActive(self, top):

        prtop = pathFolder.createFolder(self.prout + "top" + str(top) + "/")

        pafftop = prtop + "Aff" + str(top) + ".csv"
        if not path.exists(pafftop):
            # create top descriptor
            daff = toolbox.loadMatrixToDict(self.paff)
            if len(daff.keys()) <= top:
                copyfile(self.paff, pafftop)
            else:
                laff = []
                for chemID in daff.keys():
                    laff.append(float(daff[chemID]["Aff"]))
                minAff = sorted(laff, reverse=True)[top-1]
                lchem = []
                for chemID in daff.keys():
                    if float(daff[chemID]["Aff"]) >= minAff:
                        lchem.append(chemID)

                filout = open(pafftop, "w")
                filout.write("CHEMBLID\tAff\tType\n")
                for chem in lchem:
                    filout.write("%s\t%s\t%s\n"%(daff[chem]["CHEMBLID"], daff[chem]["Aff"], daff[chem]["Type"]))
                filout.close()


        runExternalSoft.clusterize(self.pdesc, pafftop, self.typeAff, self.cutoff, prtop)


    def activityCliff(self):


        runExternalSoft.activityCliff(self.pdesc, self.paff, self.typeAff, self.cutoff, self.prout)