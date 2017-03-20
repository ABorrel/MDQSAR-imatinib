from nams import nams





class MCSMatrix:

    def __init__(self, lsmiles):

        self.lsmi = lsmiles

    def computeMatrixMCS(self, pfilout):

        filout = open(pfilout, "w")

        i = 0
        imax = len(self.lsmi)
        dout = {}
        while i < imax:
            i += 1



        return


def get_Tanimoto(smile1, smile2):
    ms = nams.Nams()

    #print(smile1, smile2)
    mol_t1 = ("smi", smile1)
    mol_t2 = ("smi", smile2)
    mol1, mol_info1 = ms.get_mol_info(mol_t1[0], mol_t1[1])
    mol2, mol_info2 = ms.get_mol_info(mol_t2[0], mol_t2[1])

    #print mol_info1.keys()
    #print mol_info2.keys()
    # similarity combination
    sim12, d_atoms12 = ms.get_similarity(mol_info1, mol_info2)
    sim21, d_atoms21 = ms.get_similarity(mol_info2, mol_info1)
    sim11, d_atoms11 = ms.get_similarity(mol_info1, mol_info1)
    sim22, d_atoms22 = ms.get_similarity(mol_info2, mol_info2)

    #print sim11, d_atoms11

    #test similarity
    sizeMCS = len(d_atoms12.keys())
    nbatomdiff = max([len(d_atoms11.keys()), len(d_atoms22.keys())]) - len(d_atoms12.keys())
    text = str(sizeMCS) + "-" + str(nbatomdiff)
    #print nbatomdiff

    ##### TEST FOR NAMS SCORE #######
    #################################

    #print len(d_atoms12.keys()), len(d_atoms11.keys()), len(d_atoms22.keys())

    #ks=d_atoms12.keys()
    #ks.sort()

    #Print the atomic alignment between the molecules and its similarity score
    #for k in ks:
    #    print "\t%5d (%3s)  -%3d  (%3s) --> %6.2f" % (k[0], mol1.atoms[k[0]-1].OBAtom.GetType(),
    #                                                  k[1], mol2.atoms[k[1]-1].OBAtom.GetType(),
    #                                                  d_atoms12[k])


    # sim11
    #ks=d_atoms11.keys()
    #ks.sort()

    #print("###############################")
    #Print the atomic alignment between the molecules and its similarity score
    #for k in ks:
    #    print "\t%5d (%3s)  -%3d  (%3s) --> %6.2f" % (k[0], mol1.atoms[k[0]-1].OBAtom.GetType(),
    #                                                  k[1], mol2.atoms[k[1]-1].OBAtom.GetType(),
    #                                                  d_atoms11[k])



    # based on a Jaccard score
    #print sim12, sim21, sim11, sim22
    #print d_atoms12
    score12 = sim12 / (sim11 + sim22 - sim12)
    score21 = sim21 / (sim11 + sim22 - sim21)
    #print score12

    return [score12, text]


#a = "/home/borrel/Yue_project/result/Pi_LSR_LGDsimilarity/KS-5_1B38/cycle-CON_1N3_3ULI.smi"
#b = "/home/borrel/Yue_project/result/Pi_LSR_LGDsimilarity/KS-5_1B38/onlyC_PM1_1PYE.smi"

#get_Tanimoto(a,b)