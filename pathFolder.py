from os import listdir, remove, makedirs

PR_REF = "/home/aborrel/imitanib/"
PR_RESULT = "/home/aborrel/imitanib/results/"
PR_TEMP3D = "/home/aborrel/imitanib/results/temp3D/"
PR_COMPOUNDS = "/home/aborrel/imitanib/results/compounds/"
PR_ANALYSIS = "/home/aborrel/imitanib/results/analysis/"


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