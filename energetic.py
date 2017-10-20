import toolBox
import runOtherProg
import parseNACCESS

from os import path, system
import string


# hydrophobic -> kyte doolittle
kyte_index = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
                         "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7,
                         "W": -0.9, "Y": -1.3, "V": 4.2}

res_pos = "HKR"
res_neg = "DE"



def hydrophobicityKyte(dres, proportion=1):

    dout = {}
    dout["kyteScore"] = 0
    for resname in dres.keys():
        res = resname.split("_")[0]
        rescode = toolBox.transformAA(res)
        dout["kyteScore"] += kyte_index[rescode]

    if proportion == 1:
        dout["p_kypeScore"] = dout["kyteScore"] / len(dres.keys())

    return dout



def chargeRes(dres, proportion=1):

    dout = {}
    dout["charges"] = 0
    for resname in dres.keys():
        res = resname.split("_")[0]
        coderes = toolBox.transformAA(res)
        if coderes in res_neg:
            dout["charges"] += -1
        if coderes in res_pos:
            dout["charges"] += 1

    if proportion == 1:
        dout["p_charges"] = dout["charges"] / len(dres.keys())

    return dout


def ASADesc(dres, latoms, pprotein, ppocketatom, ppocketres):

    if not path.exists(pprotein.split(".")[0] + ".asa"):
        #compute NACCESS
        lpSA = runOtherProg.runNACESS(pprotein)
    else:
        lpSA = [pprotein.split(".")[0] + ".asa", pprotein.split(".")[0] + ".rsa"]

    print lpSA

    # extract pocket asa and rsa
    if not path.exists(ppocketatom.split(".")[0] + "_ACC.asa"):
        lpSApocket = get_ASAtom(pASA=lpSA[0], lpockatoms=latoms, ppocketatom=ppocketatom)
    else:
        lpSApocket = [ppocketatom[0:-4] + "_ACC.pdb", ppocketatom[0:-4] + "_ACC.asa"]

    #if not path.exists(ppocketres.split(".")[0] + "_ACC.rsa"):
    #    pRSApocketres = get_ASRes(pRSA=lpSA[-1], dresatom=dres, ppocketatom=ppocketatom)
    #else:
    #    pRSApocketres = ppocketres.split(".")[0] + "_ACC.pdb"

    descSA = SAScores(pSAPocket=lpSApocket[0])
    descSA["sumASA"] = accessibilityPocketbyAtom(lpSApocket[-1])

    return descSA



def get_ASRes(pRSA, ppocketatom, dresatom, SA_limit=20.0, debug=0):
    """
    Get residues exposed, exposed -> relative value
    args: -> path file rsa
          -> path file pocket
          -> value relative
    return: NONE (write file with ACC extension)
    # rm files
    os.system ("rm " + path_filin_pocket[0:-4] + "_ACC.rsa")
    os.system ("rm " + path_filin_pocket[0:-4] + "_ACC.pdb")
    """
    if path.exists(ppocketatom[0:-4] + "_ACC.rsa") and path.exists(ppocketatom[0:-4] + "_ACC.pdb"):
        return ppocketatom[0:-4] + "_ACC.pdb"

    # file out
    filoutRSA = open(ppocketatom[0:-4] + "_ACC.rsa", "w")
    filoutPDBACC = open(ppocketatom[0:-4] + "_ACC.pdb", "w")

    lres_rsa_parsed = parseNACCESS.fileRSA(pRSA)

    for resname in dresatom.keys():
        res_NUM = resname.split("_")[1]
        chain = resname.split("_")[-1]
        if debug:
            print res_NUM, "Residues num"
            print chain, "Chain pocket"
        for res_rsa in lres_rsa_parsed:
            if debug:
                print res_rsa["resSeq"], res_rsa["chainID"]
            if res_rsa["resSeq"] == res_NUM and res_rsa["chainID"] == chain:
                if res_rsa["REL"] >= SA_limit:
                    for atomres in dresatom[resname]:
                        atomres.writeAtom(filoutPDBACC)
                    filoutRSA.write(res_rsa["line"])
                    break

    filoutRSA.close()
    filoutPDBACC.close()

    return ppocketatom[0:-4] + "_ACC.pdb"


def get_ASAtom(pASA, ppocketatom, lpockatoms, SA_limit=0.0, debug=1):
    """
    Get residues exposed, exposed -> relative value
    args: -> path file rsa
          -> path file pocket
          -> value relative
    return: NONE (write file with ACC extension)
    os.system ("rm " + path_atom_pocket[0:-4] + "_ACC.asa")
    os.system ("rm " + path_atom_pocket[0:-4] + "_ACC.pdb")
    """

    filoutASA = open(ppocketatom[0:-4] + "_ACC.asa", "w")
    filoutPDBASA = open(ppocketatom[0:-4] + "_ACC.pdb", "w")

    lrsa_parsed = parseNACCESS.fileASA(pASA)

    #print lrsa_parsed
    #print dir(lpockatoms[0].chainID)
    for pockatoms in lpockatoms:
        #print pockatoms
        atom_NUM = int(pockatoms.serial)
        chain = pockatoms.chainID
        if debug:
            print chain, "Chain pocket"
            print atom_NUM, "atom num"
        for atom_asa in lrsa_parsed:
            if atom_asa["atomSeq"] == atom_NUM and chain == atom_asa["chainID"]:
                if atom_asa["ABS"] >= SA_limit:
                    if debug:
                        print "Accessible Atom"
                    pockatoms.writeAtom(filoutPDBASA)
                    filoutASA.write(atom_asa["line"])
                    break

    filoutASA.close()
    filoutPDBASA.close()
    return [ppocketatom[0:-4] + "_ACC.pdb", ppocketatom[0:-4] + "_ACC.asa"]


def SAScores (pSAPocket):

    dout = {}

    #convert mol2
    pmol2 = runOtherProg.babelPDBtoMOL2(pSAPocket)

    system("grep \" O\.\" " + pmol2 + " | wc -l > temp_O.out")
    fileO = open("temp_O.out", 'r')
    nO = fileO.readline()
    nO = int(string.replace(nO, '\n', ''))
    fileO.close()
    system("rm temp_O.out")
    system("grep \" N\.\" " + pmol2 + " | wc -l > temp_N.out")
    fileN = open("temp_N.out", 'r')
    nN = fileN.readline()
    nN = int(string.replace(nN, '\n', ''))
    fileN.close()
    system("rm temp_N.out")
    system("grep \" S\.\" " + pmol2 + " | wc -l > temp_S_all.out")

    # Without S methionine ... cf Olivier
    system("grep \" S\..*CYS\" " + pmol2 + " | wc -l > temp_S.out")
    fileS = open("temp_S.out", 'r')
    nS = fileS.readline()
    nS = int(string.replace(nS, '\n', ''))
    fileS.close()
    system("rm temp_S.out")
    fileSall = open("temp_S_all.out", 'r')
    nSall = fileSall.readline()
    nSall = int(string.replace(nSall, '\n', ''))
    fileSall.close()
    system("rm temp_S_all.out")
    system("grep \" C\.3\" " + pmol2 + " | wc -l > temp_C.out")
    fileC = open("temp_C.out", 'r')
    nC = fileC.readline()
    nC = int(string.replace(nC, '\n', ''))
    fileC.close()
    system("rm temp_C.out")

    # polarity ratio
    try:
        polarity = (nO + nN + nS) / float(nO + nN + nSall + nC)
        polarity = "%.2f" % polarity
    except:
        polarity = "NA"

    # hydrophobicity ratio voir Burgoyne et Jackson Bioinformatics
    try:
        hydroSA = (nC + nS) / float(nO + nN + nSall + nC)
        hydroSA = "%.2f" % hydroSA
    except:
        hydroSA = "NA"


    dout["polaritySA"] = polarity
    dout["hydroSA"] = hydroSA

    return dout


def accessibilityPocketbyAtom(ppockASA):
    """Retrieve sum of atom ABS
    arg: pah file asa
    return: accessibility (float)"""

    list_atom_asa = parseNACCESS.fileASA(ppockASA)
    sumASA = 0.0

    for atom_asa in list_atom_asa:
        sumASA = sumASA + atom_asa["ABS"]

    return sumASA

