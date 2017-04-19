from os import listdir, remove, makedirs

PR_REF = "/data/aborrel/imatinib/"
PR_RESULT = "/data/aborrel/imatinib/results/"
PR_TEMP3D = "/data/aborrel/imatinib/results/temp3D/"
PR_COMPOUNDS = "/data/aborrel/imatinib/results/compounds/"
PR_ANALYSIS = "/data/aborrel/imatinib/results/analysis/"


def cleanFolder(prin=PR_TEMP3D):
    lfiles = listdir(prin)
    if len(lfiles) != 0:
        for filin in lfiles:
            # problem with folder
            remove(prin + filin)
    return prin


def analyses(psub):

    if psub == "":
        return PR_ANALYSIS
    else:
        try: makedirs(PR_ANALYSIS + psub + "/")
        except: pass

    return PR_ANALYSIS + psub + "/"
