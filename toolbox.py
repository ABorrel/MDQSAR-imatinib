from re import search, sub
from shutil import copyfile
from os import listdir
from time import sleep

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


def writeFilesParamaterDesmond(pmsjout, pcfgout, timeMD, intervalFrame):

    # CFG file
    copyfile("./Desmond.cfg", pcfgout)

    filincfg = open(pcfgout, "r")
    llinescfg = filincfg.readlines()
    filincfg.close()

    # modifie CFG file
    filoutcfg = open(pcfgout, "w")
    flagtime = 0
    for linecfg in llinescfg:
        if search("^time = ", linecfg):
            filoutcfg.write("time = " + str(timeMD) + "\n")
            flagtime = 1
        elif search("^   interval = ", linecfg) and flagtime == 1:
            filoutcfg.write("   interval = " + str(intervalFrame) + "\n")
        else:
            filoutcfg.write(linecfg)
    filoutcfg.write("\n")
    filoutcfg.close()


    # MSJ file
    copyfile("./Desmond.msj", pmsjout)

    filinmsj = open(pmsjout, "r")
    llinesmsj = filinmsj.readlines()
    filinmsj.close()

    # modifie MSJ file
    filoutmsj = open(pmsjout, "w")
    for linemsj in llinesmsj:
        if search("^   cfg_file = ", linemsj):
            # need only the file name
            filoutmsj.write("   cfg_file = \"" + str(pcfgout.split("/")[-1]) + "\"\n")
        else:
            filoutmsj.write(linemsj)
    filoutmsj.write("\n")
    filoutmsj.close()


def parallelLaunch(lprout, maxJob, regexout, wait=6):

    # control number of run
    flag = 0
    while len(lprout) >= maxJob:
        i = 0
        imax = len(lprout)
        while i < imax:
            lfilesMD = listdir(lprout[i])
            for fileMD in lfilesMD:
                if search(regexout, fileMD):
                    del lprout[i]
                    imax = imax - 1
                    flag = 1
                    break
            if flag == 1:
                flag = 0
                continue
            else:
                i += 1

        if len(lprout) >= maxJob:
            print "Break " + str(wait) + "s"
            sleep(wait)

    return lprout


def loadMatrixTMalign(pmatrix):

    dmatrix = {}
    filin = open(pmatrix, "r")
    list_lines = filin.readlines()
    filin.close()

    m = 1
    for line_file in list_lines[2:5]:
        line_format = sub("[ ]{2,}", " ", line_file.strip())
        line_format = line_format.split(" ")
        dmatrix["t" + str(m)] = float(line_format[1])
        dmatrix["u" + str(m) + "1"] = float(line_format[2])
        dmatrix["u" + str(m) + "2"] = float(line_format[3])
        dmatrix["u" + str(m) + "3"] = float(line_format[4])
        m += 1

    return dmatrix



def loadTableFPI(pfileFPI):

    dout = {}
    filin = open(pfileFPI, "r")
    llineFPI = filin.readlines()
    for lineFPI in llineFPI[1:]:
        lelem = lineFPI.strip().split("\t")
        resID = lelem[0]
        lrespocket = lelem[1].split("-")
        lFPI = lelem[2].split("-")
        i = 0
        while i < len(lrespocket):
            if not resID in dout.keys():
                dout[resID] = {}
            dout[resID][lrespocket[i]] = lFPI[i]
            i += 1
    filin.close()
    return dout

def selectMinimalEnergyLigPrep(psdfin, psdfout):

    # case of only one
    filin = open(psdfin, "r")
    readfile = filin.read()
    filin.close()

    lsdf = readfile.split("$$$$\n")[:-1]


    if len(lsdf) == 1:
        copy(psdfin, psdfout)

    else:
        #find with the lower energy
        lenergy = []
        for sdfin in lsdf:
            energy = sdfin.split("> <r_lp_Energy>\n")[-1].split("\n")[0]
            print energy
            lenergy.append(float(energy))

        # take minimal energy
        ibest = lenergy.index(min(lenergy))
        print ibest
        filout = open(psdfout, "w")
        filout.write(lsdf[ibest] + "$$$$\n")
        filout.close()

    return psdfout


def parseKrakenX(pfilin):

    dout = {}
    filin = open(pfilin, "r")
    lines = filin.readlines()
    filin.close()
    lheader = lines[0].replace("  ", " ").strip().split(" ")

    i = 1
    nbcpd = len(lines)

    while i < nbcpd:
        lvalue = lines[i].replace("  ", " ").strip().split(" ")
        name = lvalue[0].split(".")[0]
        print name, i
        lvalue = lvalue[1:]

        dout[name] = {}

        if len(lvalue) != len(lheader):
            print len(lvalue), len(lheader), "Difference length"
            fff
        else:
            for j in range(0, len(lheader)):
                dout[name][lheader[j]] = lvalue[j]

        i += 1

    return dout


def loadTable(pfilin, din={}):

    filin = open(pfilin, "r")
    linesDesc = filin.readlines()
    filin.close()
    ldesc = linesDesc[0].strip().split("\t")[1:]

    nbcpd = len(linesDesc)
    i = 1
    while i < nbcpd:
        lval = linesDesc[i].strip().split("\t")
        nameCpd = lval[0]
        if not nameCpd in din:
            din[nameCpd] = {}

        j = 0
        while j < len(ldesc):
            val = lval[j+1]
            if val == "nan" or val == "inf":
                val = "NA"
            din[nameCpd][ldesc[j]] = val
            j +=1

        i+=1

    return din

def writeTableDesc(ddesc, pfilout):

    ldesc = []
    for cpdID in ddesc.keys():
        if ldesc == []:
            ldesc = ddesc[cpdID].keys()
        else:
            for desc in ddesc[cpdID].keys():
                if not desc in ldesc:
                    ldesc.append(desc)

    filout = open(pfilout, "w")
    filout.write("ID\t" + "\t".join(ldesc) + "\n")
    for cpdID in ddesc.keys():
        filout.write(str(cpdID))
        for desc in ldesc:
            try: filout.write("\t" + str(ddesc[cpdID][desc]))
            except: filout.write("\tNA")
        filout.write("\n")

    filout.close()





