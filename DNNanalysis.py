from os import listdir
from re import search

import pathFolder
import runExternalSoft




def computePerformance(prin):


    for prmodel in listdir(prin):
        dfile = {}
        print prmodel
        lfilpred = listdir(prin + prmodel + "/")
        for filepred in lfilpred:
            if search("trainpredictions.csv", filepred):
                dfile["train"] = prin + prmodel + "/" + filepred
            elif search("trainfit_testpredictions.csv", filepred):
                dfile["test"] = prin + prmodel + "/" + filepred
            elif search("CV_predictions_updated.csv", filepred):
                dfile["CV"] = prin + prmodel + "/" + filepred
        dfile["prout"] = prin + prmodel + "/"

        # compute performance
        runExternalSoft.computeRegPerf(dfile["train"], "Train", dfile["prout"])
        runExternalSoft.computeRegPerf(dfile["CV"], "CV", dfile["prout"])
        runExternalSoft.computeRegPerf(dfile["test"], "Test", dfile["prout"])



###################
#  Deep learning  #
###################

PR_DL = pathFolder.PR_HOME + "imatinib/DL/Split/"


computePerformance(PR_DL)
