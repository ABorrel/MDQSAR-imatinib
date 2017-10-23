import composition
import PDB
import energetic
import geometric

from os import path

class pocket:

    def __init__(self, ppocket, pPDB):

        self.pprotein = pPDB
        self.ppocket = ppocket

        cPocket = PDB.PDB(ppocket)  # not included hydrogen
        self.latoms = cPocket.get_lAtoms()
        self.byresall = cPocket.get_byres(onlyres=1)



    def get_alldescs(self):
        self.get_compo(1)
        self.get_energetic(1)
        self.get_geometric()

        self.ldescs = self.lcompo + self.lenergy + self.lgeo


    def get_compo(self, proportion=1):

        dcompo = {}
        dcompo.update(composition.compoRes(self.byresall, proportion))
        dcompo.update(composition.compoResType(self.byresall, proportion))
        dcompo.update(composition.compoAtom(self.latoms, proportion))
        self.compo = dcompo
        self.lcompo = dcompo.keys()

    def get_energetic(self, proportion=1):

        denergy = {}
        denergy.update(energetic.hydrophobicityKyte(self.byresall, proportion))
        denergy.update(energetic.chargeRes(self.byresall, proportion))
        denergy.update(energetic.ASADesc(self.byresall, self.latoms, self.pprotein, self.ppocket, self.byresall))
        self.energy = denergy
        self.lenergy = denergy.keys()


    def get_geometric(self):

        dgeo = {}
        dgeo.update(geometric.get_RADII(self.ppocket))
        self.geometry = dgeo
        self.lgeo = dgeo.keys()



    def writeDesc(self, pdesc, rowname):

        if not path.exists(pdesc):
            fildesc = open(pdesc, "w")
            #header
            fildesc.write("\t".join(self.ldescs) + "\n")
        else:
            fildesc = open(pdesc, "a")

        fildesc.write(str(rowname))

        for desc in self.ldescs:
            if desc in self.compo.keys():
                val = self.compo[desc]
            elif desc in self.energy.keys():
                val = self.energy[desc]
            elif desc in self.geometry.keys():
                val = self.geometry[desc]
            else: print "ERROR"
            fildesc.write("\t" + str(val))
        fildesc.write("\n")
        fildesc.close()


