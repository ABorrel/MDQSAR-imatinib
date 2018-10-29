from os import path, listdir

import ligand
import runExternalSoft
import pathFolder

class DescriptorLig:
    def __init__(self, lccpd, prPoses, corVal, maxquantile, compute1D2D, compute3D, prout):

        self.lccpd = lccpd
        self.prPoses = prPoses
        self.compute1D2D = compute1D2D
        self.compute3D = compute3D
        self.corVal = corVal
        self.maxquantile = maxquantile
        self.prout = prout


    def computeDesc(self):

        self.pdesc = {}

        if self.compute1D2D == 1:
            pdesc1D2D = self.prout + "1D2D.csv"
            self.pdesc["1D2D"] = pdesc1D2D
            if path.exists(pdesc1D2D) and path.getsize(pdesc1D2D) > 10:
                compute1D2D = 0
            else:
                compute1D2D = 1
        else:
            compute1D2D = 0

        if self.compute3D == 1:
            pdesc3D = self.prout + "3D.csv"
            self.pdesc["3D"] = pdesc3D
            if path.exists(pdesc3D) and path.getsize(pdesc3D) > 10:
                compute3D = 0
            else:
                compute3D = 1
        else:
            compute3D = 0

        if compute1D2D == 1 and compute3D == 1:
            p1D2D3D = self.prout + "1D2D3D.csv"
            self.pdesc["1D2D3D"] = p1D2D3D

        if compute1D2D == 1 or compute3D == 1:
            logfile = open(self.prout + "desc.log", "w")
            for compound in self.lccpd:
                cDesc = ligand.Descriptors(compound, logfile)
                if compute3D == 1:
                    CHEMBLID = cDesc.compound[cDesc.kID]
                    psdf3D = self.prPoses + CHEMBLID + ".1.sdf"
                    if not path.exists(psdf3D):
                        continue

                if compute1D2D == 1:
                    cDesc.get_descriptor1D2D()
                    cDesc.writeTablesDesc(pdesc1D2D, "1D2D")
                if compute3D == 1:
                    cDesc.get_descriptor3D(self.prPoses)
                    cDesc.writeTablesDesc(pdesc3D, "3D")
                if compute3D == 1 and compute1D2D == 1:
                    cDesc.writeTablesDesc(p1D2D3D, "1D2D3D")


    def dendoAffinity(self, typeDesc, paff):

        prout = pathFolder.createFolder(self.prout + "dendoAff_" + str(typeDesc) + "/")
        pdesc = self.pdesc[typeDesc]
        if not path.exists(pdesc):
            print "No path descriptor", pdesc, "found"
            return 1

        runExternalSoft.DescAnalysis(pdesc, paff, prout, self.corVal, self.maxquantile, PCA=0, corMatrix=0, hist=0, dendo=1, clustering=0)

        return 0


def generatePNG():

    lsmi = listdir(pathFolder.PR_SMI)
    for smi in lsmi:
        psmi = pathFolder.PR_SMI + smi
        ppng = pathFolder.PR_PNG + smi[:-3] + "png"
        runExternalSoft.molconvert(psmi, ppng)
