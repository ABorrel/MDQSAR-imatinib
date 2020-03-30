from numpy import mean

import toolbox


class RMSD:
    def __init__(self, prin):
        self.prin = prin

    def loadRMSDs(self, lload = ["ligand", "protein", "residues"]):

        for typeRMSD in lload:
            if typeRMSD == "ligand":
                dlig = {}
                dlig["RMSF"] = {}

                pRMSF = self.prin + typeRMSD + "/ligRMSF"
                dlig["RMSF"] = toolbox.loadMatrixToDict(pRMSF)

                pShaEP = self.prin + typeRMSD + "/ligShaEP"
                dlig["ShaEP"] = toolbox.matrixToList(pShaEP)

                self.lig = dlig

            if typeRMSD == "protein":
                pRMSD = self.prin + typeRMSD + "/protRMSD"
                dprot = toolbox.matrixToList(pRMSD)
                self.prot = dprot

            if typeRMSD == "residues":
                pRes = self.prin + typeRMSD + "/resRMSD"
                dres = toolbox.loadMatrixToDict(pRes)
                self.res = dres


    def MRMSDprot(self):
        if not "prot" in self.__dict__:
            self.loadRMSDs(["protein"])

        dout = {}
        for frame in self.prot:
            for k in frame.keys():
                if k != "Time":
                    if not k in dout.keys():
                        dout[k] = []

                    dout[k].append(float(frame[k]))

        return [mean(dout['RMSDC']), mean(dout['RMSDall']), mean(dout['Dmax'])]


    def MRMSDlig(self):

        if not "lig" in self.__dict__:
            self.loadRMSDs(["ligand"])
        lRMSD = []
        for atom in self.lig["RMSF"].keys():
            lRMSD.append(float(self.lig["RMSF"][atom]["RMSF"]))

        return mean(lRMSD)





