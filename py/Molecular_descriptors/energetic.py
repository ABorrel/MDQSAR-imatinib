import toolbox
import runExternalSoft
import PDB
import parseRSAfile
#import parseNACCESS

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
        rescode = toolbox.transformAA(res)
        dout["kyteScore"] += kyte_index[rescode]

    if proportion == 1:
        dout["p_kypeScore"] = dout["kyteScore"] / len(dres.keys())

    return dout



def chargeRes(dres, proportion=1):

    dout = {}
    dout["charges"] = 0
    for resname in dres.keys():
        res = resname.split("_")[0]
        coderes = toolbox.transformAA(res)
        if coderes in res_neg:
            dout["charges"] += -1
        if coderes in res_pos:
            dout["charges"] += 1

    if proportion == 1:
        dout["p_charges"] = dout["charges"] / len(dres.keys())

    return dout


def ASADesc(pprotein, ppocketatom):

    dout = {}
    # run freeSASA for the protein
    lfilout = runExternalSoft.runFreeSASA(pprotein, pprotein[0:-4] + "_asa.pdb", rsa=1)

    dASA = ASAHydrophobicityPolarity(lfilout[0], ppocketatom)
    dout["polarityASA"] = dASA[0]
    dout["hydrophobicityASA"] = dASA[1]
    dout.update(sumAccesbilityBS(ppocketatom, lfilout[0], lfilout[1]))

    return dout



def ASAHydrophobicityPolarity(ppdbasa, pBS):

    cPDBasa = PDB.PDB(ppdbasa)
    latomASA = cPDBasa.get_lAtoms()
    cBS = PDB.PDB(pBS)
    latomBS = cBS.get_lAtoms()

    dcompute = {}
    dcompute["C"] = []
    dcompute["O"] = []
    dcompute["N"] = []
    dcompute["Scys"] = []
    dcompute["Smet"] = []

    for atomBS in latomBS:
        for atomProt in latomASA:
            if atomBS.chainID == atomProt.chainID and atomBS.name == atomProt.name and atomBS.resName == atomProt.resName and atomBS.serial == atomProt.serial:
                if atomBS.element != "S":
                    dcompute[atomBS.element].append(float(atomProt.Bfact))
                else:
                    if atomBS.resName == "CYS":
                        dcompute["Scys"].append(float(atomProt.Bfact))
                    elif atomBS.resName == "MET":
                        dcompute["Smet"].append(float(atomProt.Bfact))
                    else:
                        print atomBS
                        dddd
                break

    polarityASA = (sum(dcompute["O"]) + sum(dcompute["N"]) + sum(dcompute["Scys"])) / (sum(dcompute["O"]) + sum(dcompute["N"]) + sum(dcompute["Scys"]) + sum(dcompute["Smet"]) + sum(dcompute["C"]))
    hydrophobicityASA = (sum(dcompute["C"]) + sum(dcompute["Smet"])) / (sum(dcompute["O"]) + sum(dcompute["N"]) + sum(dcompute["Scys"]) + sum(dcompute["Smet"]) + sum(dcompute["C"]))

    return [polarityASA, hydrophobicityASA]



def sumAccesbilityBS(ppocketatom, proteinASA, proteinRSA):


    cASA = PDB.PDB(proteinASA)
    cASA.get_lAtoms()

    cRSA = parseRSAfile.RSA(proteinRSA)

    cBS = PDB.PDB(ppocketatom)
    latomBS = cBS.get_lAtoms()

    dout = {}
    dout["sumASA"] = 0.0
    for atomBS in latomBS:
        for atomProt in cASA.latom:
            if atomBS.chainID == atomProt.chainID and atomBS.name == atomProt.name and atomBS.resName == atomProt.resName and atomBS.serial == atomProt.serial:
                dout["sumASA"] += float(atomProt.Bfact)
                break

    dout["sumRSAabs"] = 0.0
    dout["sumRSArel"] = 0.0
    for res in cRSA.lres:
        for atomBS in latomBS:
            if atomBS.chainID == res.chainID and atomBS.resName == res.resName and atomBS.resSeq == res.resSeq:
                if res.ABSall != "N/A":
                    dout["sumRSAabs"] += float(res.ABSall)
                if res.RELall != "N/A":
                    dout["sumRSArel"] += float(res.RELall)
                break
    return dout