from os import listdir, remove, makedirs, path
from shutil import rmtree


# for personal computer
#PR_REF = "/home/aborrel/imitanib/"
#PR_RESULT = "/home/aborrel/imitanib/results/"
#PR_TEMP3D = "/home/aborrel/imitanib/results/temp3D/"
#PR_COMPOUNDS = "/home/aborrel/imitanib/results/compounds/"
#PR_ANALYSIS = "/home/aborrel/imitanib/results/analysis/"

# persomal NIH
PR_REF = "/home/borrela2/imatinib/"
PR_RESULT = "/home/borrela2/imatinib/results/"
PR_TEMP3D = "/home/borrela2/imatinib/results/temp3D/"
PR_COMPOUNDS = "/home/borrela2/imatinib/results/compounds/"
PR_ANALYSIS = "/home/borrela2/imatinib/results/analysis/"

# for monster
#PR_REF = "/data/aborrel/imatinib/"
#PR_RESULT = "/data/aborrel/imatinib/results/"
#PR_TEMP3D = "/data/aborrel/imatinib/results/temp3D/"
#PR_COMPOUNDS = "/data/aborrel/imatinib/results/compounds/"
#PR_ANALYSIS = "/data/aborrel/imatinib/results/analysis/"


def cleanFolder(prin=PR_TEMP3D):
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


def analyses(psub):

    if psub == "":
        return PR_ANALYSIS
    else:
        try: makedirs(PR_ANALYSIS + psub + "/")
        except: pass

    return PR_ANALYSIS + psub + "/"
