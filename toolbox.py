from re import search


def transformAA(aa):
    aa = aa.upper()
    dico_code = {"S": "SER", "T": "THR", "N": "ASN", "Q": "GLN", "E": "GLU", "D": "ASP", "K": "LYS", "R": "ARG",
                 "H": "HIS", "M": "MET", "C": "CYS", "W": "TRP", "F": "PHE", "Y": "TYR", "A": "ALA", "V": "VAL",
                 "L": "LEU", "I": "ILE", "P": "PRO", "G": "GLY"}

    if len(aa) == 1:
        return dico_code[aa]
    else:
        for aa_one in dico_code.keys():
            if dico_code[aa_one] == aa:
                return aa_one


def checkENDFinalLinePDBfile(pinPDB):
    """
    Check if character END in end file PDB
    args: -> path file PDB
    return : -> NONE change directly file PDB
    """
    filin = open(pinPDB, "r")
    list_lines = filin.readlines()
    filin.close()
    if not search("^END", list_lines[-1]):
        if (list_lines[-1][-1] == "\n"):
            filout = open(pinPDB, "a")
            filout.write("END\n")
            filout.close()
        else:
            filout = open(pinPDB, "a")
            filout.write("\nEND\n")
            filout.close()
    elif not search("^END\n", list_lines[-1]):
        filout.write("\n")
        filout.close()
    else:
        pass


def checkHEADERinitialLinePDBfile(pinPDB):
    """
    Check if character HEADER in beginning file PDB
    args: -> path file PDB
    return : -> NONE change directly file PDB
    """

    filin = open(pinPDB, "r")
    element = filin.read()
    filin.close()
    if not search("^HEADER", element):
        filout = open(pinPDB, "w")
        filout.write("HEADER \n")
        filout.write(element)
        filout.close()





def mergeDict(l_dict):

    dout = {}

    lk = l_dict[0].keys()

    for k in lk:
        lval = [l_dict[i][k] for i in range(0,len(l_dict))]
        #identic
        if lval.count(lval[0]) == len(lval):
            dout[k] = lval[0]
        else:
            dout[k] = "----".join(lval)

    return dout




def convertUnit(l_values, l_units):
    """Convert list of affinity un uM"""

    #print l_values
    #print l_units
    lout = []
    i = 0
    imax = len(l_values)
    while i < imax:
        if l_units[i] == "uM" or l_units[i] == "10'-6M" or l_units[i] == "umol/L" or l_units[i] == "microM":
            lout.append(l_values[i])
            i += 1
        elif l_units[i] == "nM":
            val = float(l_values[i])
            val = val/1000
            lout.append(val)
            i += 1
        elif l_units[i] == "10'-3microM":
            val = float(l_values[i])
            val = val/1000
            lout.append(val)
            i += 1

        else:
            print "sssss",l_units[i]
            ffff

    #print lout

    return lout




def loadMatrix(pfilin):
    "Case of square matrix"

    filin = open(pfilin, "r")
    llinesMatrix = filin.readlines()
    filin.close()

    dout = {}
    lcompID = llinesMatrix[0].strip().split("\t")
    nbComp = len(lcompID)

    i = 1
    while i < nbComp:
        if not lcompID[i] in dout.keys():
            dout[lcompID[i]] = {}
        j = 1
        lval = llinesMatrix[i].strip().split("\t")
        while j < nbComp:
            dout[lcompID[i]][lcompID[j]] = lval[j]
            j += 1
        i += 1

    return dout


