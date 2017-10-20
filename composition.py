import string
from numpy import mean

import toolBox

# Residues -> VENN Faure 2008
res_all = "ACDEFGHIKLMNPQRSTVWY"
res_pos = "HKR"
res_neg = "DE"
res_aromatic = "FYHW"
res_polar = "CDEHKNQRSTWY"
res_tiny = "ACGS"
res_hydrophobe = "GACTVLIMFWYHK"
res_aliphatic = "ILV"
res_small = "CVTGASDNP"
res_charged = "DERKH"

# res_Hacceptors = "DEHNQ"
# res_Hdonors = "CHKNQRSTWY"



# Atoms
Atom = "CNOS"

# reference atom type (Milletti 2010) -> table 1 for aromatic carbon
Ctype = {"R": ["CB", "CD", "CG"], "M": ["CB"], "F": ["CB"], "L:": ["CB", "CD1", "CD2", "CG"], "W": ["CB"],
               "D": ["CB"], "K": ["CB", "CE"], "H": ["CB"], "V": ["CB", "CG1", "CG2"], "Q": ["CB", "CG", "CD"],
               "A": ["CB"], "E": ["CB"], "P": ["CB", "CD", "CG"], "C": ["CB"], "Y": ["CB"], "N": ["CB"],
               "I": ["CB", "CD1", "CG1", "CG2"]}
Car = {"F": ["CD1", "CD2", "CE1", "CE2", "CG", "CZ"],
                 "W": ["CD1", "CD2", "CE2", "CE3", "CG", "CH2", "CZ2", "CZ3"], "H": ["CD2", "CE1"],
                 "Y": ["CD1", "CD2", "CE1", "CE2", "CG", "CZ"]}
Carg = {"S": ["CB"], "T": ["CB"]}
Ntype = {"A": ["N"], "C": ["N"], "D": ["N"], "E": ["N"], "F": ["N"], "G": ["N"], "H": ["N"], "I": ["N"],
               "K": ["N"], "L": ["N"], "M": ["N"], "N": ["N", "ND2"], "P": ["N"], "Q": ["N", "NE2"], "R": ["N"],
               "S": ["N"], "T": ["N"], "V": ["N"], "W": ["N"], "Y": ["N"]}
ND1 = {"H": ["ND1"]}
NE2 = {"H": ["NE2"]}
Nlys = {"K": ["NZ"]}
Ntrp = {"W": ["NE1"]}
Otype = {"A": ["O"], "C": ["O"], "D": ["O"], "E": ["O"], "F": ["O"], "G": ["O"], "H": ["O"], "I": ["O"],
               "K": ["O"], "L": ["O"], "M": ["O"], "N": ["O"], "P": ["O"], "Q": ["O"], "R": ["O"], "S": ["O"],
               "T": ["O"], "V": ["O"], "W": ["O"], "Y": ["O"]}
Ocoo = {"E": ["CG"]}
Ooh = {"S": ["CB"], "T": ["CB"]}
Otyr = {"Y": ["OH"]}
Stype = {"C": ["SG"]}
Ccoo = {"D": ["CB"], "E": ["CG"]}
Cgln = {"Q": ["CG"], "N": ["CB"]}
Hyd = {"R": ["CZ"], "M": ["SD"], "F": ["CG", "CZ"], "L": ["CG"], "W": ["CE3", "CG", "CZ2"], "H": ["CG"],
                 "V": ["CB"], "P": ["CG"], "C": ["SG"], "Y": ["CG", "CZ"]}

# old version
dhydrophobic = {"R": ["CZ"], "M": ["SD"], "F": ["CG", "CZ"], "L": ["CG"], "W": ["CE3", "CG", "CZ2"],
                         "H": ["CG"], "V": ["CB"], "P": ["CG"], "C": ["SG"], "Y": ["CG", "CZ"], "I": ["CB"]}
dHacceptor = {"D": ["OD1", "OD2"], "E": ["OE1", "OE2"], "H": ["CG", "CD2", "NE2", "CE1", "ND1"], "N": ["OD1"],
                  "Q": ["OE1"], "ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"], "HIS": ["CG", "CD2", "NE2", "CE1", "ND1"],
                  "ASN": ["OD1"], "GLN": ["OE1"]}
dHdonor = {"C": ["SG"], "H": ["CG", "CD2", "NE2", "CE1", "ND1"], "K": ["NZ"], "N": ["ND2"], "Q": ["NE2"],
               "R": ["NE", "NH1", "NH2"], "S": ["OG"], "T": ["OG1"], "Y": ["OH"], "W": ["NE1"]}






def compoRes(drespocket, proportion=1):

    dout = {}
    dout["CRES"] = len(drespocket.keys())
    for res in res_all:
        dout[res] = 0.0

    #count global
    for resPocket in drespocket.keys():
        nameres = resPocket.split("_")[0]
        dout[toolBox.transformAA(nameres)] += 1

    #proportion
    if proportion == 1:
        for res in res_all:
            dout["p" + res] = dout[res]/dout["CRES"]
    return dout



def compoResType(drespocket, proportion=1):

    cres = len(drespocket.keys())
    dout = {}
    ltyperes = ["pos", "neg", "arom", "polar", "tiny", "hydrophobe", "aliphatic", "small", "charged"]

    for typeres in ltyperes:
        dout[typeres] = 0.0

    for resPocket in drespocket.keys():
        nameres = resPocket.split("_")[0]
        coderes = toolBox.transformAA(nameres)

        if coderes in res_pos:
            dout["pos"] += 1
        if coderes in res_neg:
            dout["neg"] = dout["neg"] + 1
        if coderes in res_aromatic:
            dout["arom"] = dout["arom"] + 1
        if coderes in res_polar:
            dout["polar"] = dout["polar"] + 1
        if coderes in res_tiny:
            dout["tiny"] = dout["tiny"] + 1
        if coderes in res_hydrophobe:
            dout["hydrophobe"] = dout["hydrophobe"] + 1
        if coderes in res_aliphatic:
            dout["aliphatic"] = dout["aliphatic"] + 1
        if coderes in res_small:
            dout["small"] = dout["small"] + 1
        if coderes in res_charged:
            dout["charged"] = dout["charged"] + 1

    if proportion == 1:
        for typeres in ltyperes:
            dout["p" + typeres] = dout[typeres]/cres

    return dout


def compoAtom(latomspocket, proportion=1):

    dout = {}
    ltypeatom = ["Ctype", "Car", "Carg", "Ntype", "ND1", "NE2", "Nlys", "Ntrp", "Otype",  "Ocoo", "Ooh", "Otyr", "Stype",
                 "Ccoo", "Cgln", "Hyd"]
    lelem = ["C", "N", "O", "S"]

    catom = len(latomspocket)
    dout["CATOM"] = catom

    for typeatom in ltypeatom:
        dout[typeatom] = 0.0
    for elem in lelem:
        dout[elem] = 0.0

    for atompocket in latomspocket:
        if atompocket.name in Ctype:
            dout["Ctype"] = dout["Ctype"] + 1
        if atompocket.name in Car:
            dout["Car"] = dout["Car"] + 1
        if atompocket.name in Carg:
            dout["Carg"] = dout["Carg"] + 1
        if atompocket.name in Ntype:
            dout["Ntype"] = dout["Ntype"] + 1
        if atompocket.name in ND1:
            dout["ND1"] = dout["ND1"] + 1
        if atompocket.name in NE2:
            dout["NE2"] = dout["NE2"] + 1
        if atompocket.name in Nlys:
            dout["Nlys"] = dout["Nlys"] + 1
        if atompocket.name in Ntrp:
            dout["Ntrp"] = dout["Ntrp"] + 1
        if atompocket.name in Otype:
            dout["Otype"] = dout["Otype"] + 1
        if atompocket.name in Ocoo:
            dout["Ocoo"] = dout["Ocoo"] + 1
        if atompocket.name in Ooh:
            dout["Ooh"] = dout["Ooh"] + 1
        if atompocket.name in Otyr:
            dout["Otyr"] = dout["Otyr"] + 1
        if atompocket.name in Stype:
            dout["Stype"] = dout["Stype"] + 1
        if atompocket.name in Ccoo:
            dout["Ccoo"] = dout["Ccoo"] + 1
        if atompocket.name in Cgln:
            dout["Cgln"] = dout["Cgln"] + 1
        if atompocket.name in Hyd:
            dout["Hyd"] = dout["Hyd"] + 1

        if atompocket.element == "C":
            dout["C"] = dout["C"] + 1
        if atompocket.element == "O":
            dout["O"] = dout["O"] + 1
        if atompocket.element == "N":
            dout["N"] = dout["N"] + 1
        if atompocket.element == "S":
            dout["S"] = dout["S"] + 1

    if proportion == 1:
        for typeatom in ltypeatom:
            dout["p" + typeatom] = dout[typeatom] / catom

        for elem in lelem:
            dout["p" + elem] = dout[elem] / catom

    return dout






