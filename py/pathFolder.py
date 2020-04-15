from os import listdir, remove, makedirs, path
from shutil import rmtree




#PR_HOME = "C:/Users/Aborrel/research/NCSU/" #Laptop
#PR_HOME = "/home/borrela2/" #NIEHS
#PR_HOME = "/home/aborrel/" #NCSU
#PR_HOME = "/data/aborrel/" #Monster

#PR_REF = PR_HOME + "imatinib-MD/"
#PR_RESULT = PR_HOME + "imatinib/results/"
#PR_TEMP3D = PR_HOME + "imatinib/results/temp3D/"
#PR_SMI = PR_HOME + "imatinib/results/SMI/"
#PR_PNG = PR_HOME + "imatinib/results/PNG/"
#PR_ANALYSIS = PR_HOME + "imatinib/results/analysis/"


def cleanFolder(prin):
    lfiles = listdir(prin)
    if len(lfiles) != 0:
        for filin in lfiles:
            # problem with folder
            try: remove(prin + filin)
            except: rmtree(prin + filin)

    return prin


def createFolder(prin, clean=0):

    if not path.exists(prin):
        makedirs(prin)

    if clean == 1:
        cleanFolder(prin)

    return prin


#def analyses(psub):

#    if psub == "":
#        return PR_ANALYSIS
#    else:
#        try: makedirs(PR_ANALYSIS + psub + "/")
#        except: pass

#    return PR_ANALYSIS + psub + "/"




###########
# define directory for analysis
###########

PR_ROOT = "./../../"
PR_RESULT = createFolder(PR_ROOT + "results/")
PR_DATA = createFolder(PR_ROOT + "data/")
