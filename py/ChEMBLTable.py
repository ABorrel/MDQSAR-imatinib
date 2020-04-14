from copy import deepcopy
from numpy import mean, std
from os import path
from re import search
import collections
import sys

import runExternalSoft
import toolbox
import pathFolder

# import descriptor scripts
sys.path.insert(0, "./Molecular_descriptors/") # for window dev
import ligand

class ChEMBLTable:
    def __init__(self, pfilin, pr_out, ltypeAff, lBAout):
        self.pfilin = pfilin
        self.pr_out = pr_out
        self.ltypeAff = ltypeAff
        self.outlier = lBAout


    def parseCHEMBLFile(self):

        lout = []
        filin = open(self.pfilin, "r")
        llines = filin.readlines()
        filin.close()

        lhead = llines[0].split("\t")[0:-1]
        #print lhead

        i = 1
        while i < len(llines):
            lelem = llines[i].split("\t")
            #print lelem, i
            #print len(lelem), len(lhead), "***"
            #print lelem
            dout = {}
            j = 0
            while j < len(lhead):
                dout[lhead[j]] = lelem[j]
                j += 1
            lout.append(dout)
            i += 1

        self.table = lout



    def getOnlyExactConstant(self):
        """Remove inhibition in percentage and not exact value"""

        if not "table" in dir(self):
            self.parseCHEMBLFile()


        i = 0
        imax = len(self.table)
        while i < imax:
            row = self.table[i]
            if row["RELATION"] != "=":
                del self.table[i]
                imax = imax - 1
                continue
            elif row["PCHEMBL_VALUE"] == "":
                del self.table[i]
                imax = imax - 1
                continue
            elif row["PUBLISHED_UNITS"] == "":
                del self.table[i]
                imax = imax - 1
                continue
            else:
                i += 1

        #print len(self.table)

    def getByTypeOfAff(self, ltypeAff =["IC50"]):

        if not "table" in dir(self):
            self.parseCHEMBLFile()

        i = 0
        imax = len(self.table)
        while i < imax:
            row = self.table[i]
            if not row["PUBLISHED_TYPE"] in ltypeAff:
                del self.table[i]
                imax = imax - 1
                continue
            else:
                i += 1


    def MergeIdenticCHEMBLIDforACtivity(self):
        """Control quality of constant available and compute mean"""

        if not "table" in dir(self):
            self.parseCHEMBLFile()

        d_CHEMBLID = {}
        for row in self.table:
            if not row["CMPD_CHEMBLID"] in d_CHEMBLID.keys():
                d_CHEMBLID[row["CMPD_CHEMBLID"]] = []
            d_CHEMBLID[row["CMPD_CHEMBLID"]].append(deepcopy(row))

        # case no copy
        if len(d_CHEMBLID.keys()) == len(self.table):
            return

        i = 0
        imax = len(d_CHEMBLID.keys())
        while i < imax:
            #print i, d_CHEMBLID.keys()[i]
            # case not problem
            #print d_CHEMBLID.keys()[i]
            if len(d_CHEMBLID[d_CHEMBLID.keys()[i]]) == 1:
                i += 1
                continue
            else:
                # Control the published value
                l_PUBLISHED_VALUE = [float(d_CHEMBLID[d_CHEMBLID.keys()[i]][k]["PUBLISHED_VALUE"]) for k in range(0, len(d_CHEMBLID[d_CHEMBLID.keys()[i]]))]
                l_PUBLISHED_UNITS = [d_CHEMBLID[d_CHEMBLID.keys()[i]][k]["PUBLISHED_UNITS"] for k in range(0, len(d_CHEMBLID[d_CHEMBLID.keys()[i]]))]
                l_PUBLISHED_TYPE = [d_CHEMBLID[d_CHEMBLID.keys()[i]][k]["STANDARD_TYPE"] for k in range(0, len(d_CHEMBLID[d_CHEMBLID.keys()[i]]))]


                if not l_PUBLISHED_UNITS.count(l_PUBLISHED_UNITS[0]) == len(l_PUBLISHED_UNITS):
                    l_PUBLISHED_VALUE = toolbox.convertUnit(l_PUBLISHED_VALUE, l_PUBLISHED_UNITS)
                    #print l_PUBLISHED_UNITS


                # combine if same type of affinity
                # extract only most present type
                print l_PUBLISHED_TYPE
                if len(list(set(l_PUBLISHED_TYPE))) != 1:
                    if "Ki" in l_PUBLISHED_TYPE:
                        typetoextract = "Ki" # favorize most common Ki
                    else:
                        typetoextract = collections.Counter(l_PUBLISHED_TYPE).most_common()[0][0]
                    t = 0
                    tmax = len(l_PUBLISHED_TYPE)
                    while t < tmax:

                        if not l_PUBLISHED_TYPE[t] == typetoextract:
                            del l_PUBLISHED_TYPE[t]
                            del l_PUBLISHED_VALUE[t]
                            del l_PUBLISHED_UNITS[t]
                            del d_CHEMBLID[d_CHEMBLID.keys()[i]][t]
                            tmax = tmax-1
                        else:
                            t += 1



                MPUBLISHED_VALUE = mean(l_PUBLISHED_VALUE)
                SDPUBLISHED_VALUE = std(l_PUBLISHED_VALUE)

                magnitudeval = len(str(int(min(l_PUBLISHED_VALUE))))
                magnitudeSD = len(str(int(SDPUBLISHED_VALUE)))

                #print magnitudeSD, magnitudeval, "l110", d_CHEMBLID.keys()[i]

                if magnitudeval != magnitudeSD:
                    del d_CHEMBLID[d_CHEMBLID.keys()[i]]
                    imax = imax - 1
                    continue

                else:

                    l_STANDARD_VALUE = [float(d_CHEMBLID[d_CHEMBLID.keys()[i]][k]["STANDARD_VALUE"]) for k in
                                        range(0, len(d_CHEMBLID[d_CHEMBLID.keys()[i]]))]
                    l_PCHEMBL_VALUE = [float(d_CHEMBLID[d_CHEMBLID.keys()[i]][k]["PCHEMBL_VALUE"]) for k in
                                       range(0, len(d_CHEMBLID[d_CHEMBLID.keys()[i]]))]
                    MSTANDARD_VALUE = mean(l_STANDARD_VALUE)


                    d_CHEMBLID[d_CHEMBLID.keys()[i]] = [toolbox.mergeDict(d_CHEMBLID[d_CHEMBLID.keys()[i]])]
                    d_CHEMBLID[d_CHEMBLID.keys()[i]][0]["STANDARD_VALUE"] = str(MSTANDARD_VALUE)
                    d_CHEMBLID[d_CHEMBLID.keys()[i]][0]["PCHEMBL_VALUE"] = str(mean(l_PCHEMBL_VALUE))
                    d_CHEMBLID[d_CHEMBLID.keys()[i]][0]["PUBLISHED_VALUE"] = str(MPUBLISHED_VALUE)
                    i += 1

        # reformate the table
        self.table = []
        for k in d_CHEMBLID.keys():
            self.table.append(d_CHEMBLID[k][0])



    def delIdenticCHEMBLIDByCuration(self):
        """Prefered expert curation or most recent publication"""

        if not "table" in dir(self):
            self.parseCHEMBLFile()

        d_CHEMBLID = {}
        for row in self.table:
            if not row["CMPD_CHEMBLID"] in d_CHEMBLID.keys():
                d_CHEMBLID[row["CMPD_CHEMBLID"]] = []
            d_CHEMBLID[row["CMPD_CHEMBLID"]].append(deepcopy(row))

        # case no copy
        if len(d_CHEMBLID.keys()) == len(self.table):
            return

        for CHEMBLID in d_CHEMBLID.keys():
            if len(d_CHEMBLID[CHEMBLID]) == 1:
                continue
            else:
                l_curation = [d_CHEMBLID[CHEMBLID][i]["CURATED_BY"] for i in range(0,len(d_CHEMBLID[CHEMBLID]))]
                l_years = [d_CHEMBLID[CHEMBLID][i]["YEAR"] for i in range(0,len(d_CHEMBLID[CHEMBLID]))]

                # Curation
                if "Expert" in l_curation:
                    i = 0
                    imax = len(l_curation)
                    while i < len(l_curation):
                        if d_CHEMBLID[CHEMBLID][i]["CURATED_BY"] != "Expert":
                            del self.table[self.table.index(d_CHEMBLID[CHEMBLID][i])]
                            del l_years[i]
                            del l_curation[i]
                            del d_CHEMBLID[CHEMBLID][d_CHEMBLID[CHEMBLID].index(d_CHEMBLID[CHEMBLID][i])]
                            imax = imax -1
                            continue
                        i += 1

                #years
                if len(l_years) == 1:
                    continue
                else:
                    recent = max(l_years)
                    i = 0
                    imax = len(l_years)
                    while i < imax:
                        if d_CHEMBLID[CHEMBLID][i]["YEAR"] != str(recent):
                            del self.table[self.table.index(d_CHEMBLID[CHEMBLID][i])]
                            del l_years[i]
                            del d_CHEMBLID[CHEMBLID][d_CHEMBLID[CHEMBLID].index(d_CHEMBLID[CHEMBLID][i])]
                            imax -= 1
                            continue
                        i += 1

                # case where several identic - take first
                if len(d_CHEMBLID[CHEMBLID]) != 1:
                    for i in range (1,len(d_CHEMBLID[CHEMBLID])) :
                        del self.table[self.table.index(d_CHEMBLID[CHEMBLID][i])]


                #print l_years
                recent = max(l_years)

                #print recent, "eeeee"


    def selectConfidencecore(self, cutoff = 0):

        if not "table" in dir(self):
            self.parseCHEMBLFile()

        i = 0
        imax = len(self.table)
        while i < imax:
            score = float(self.table[i]["CONFIDENCE_SCORE"])
            if score < cutoff:
                del self.table[i]
                imax = imax - 1
                continue
            else:
                i += 1


    def selectAssayType(self, assaytypeSelected):

        if not "table" in dir(self):
            self.parseCHEMBLFile()

        i = 0
        imax = len(self.table)
        while i < imax:
            assaytype = self.table[i]["ASSAY_TYPE"]
            if assaytype != assaytypeSelected:
                del self.table[i]
                imax = imax - 1
                continue
            else:
                i += 1



    def removeBA(self, lBAout, debug=1):

        if not "table" in dir(self):
            self.parseCHEMBLFile()

        if debug==1: print lBAout

        i = 0
        imax = len(self.table)
        while i < imax:
            assayID = self.table[i]["ASSAY_CHEMBLID"]
            if assayID in lBAout:
                del self.table[i]
                imax = imax - 1
                continue
            else:
                i += 1




    def writeTable(self, pfilout):

        # define order to writing
        lcol = ['ACTIVITY_ID', 'CELL_ID', 'PUBLISHED_RELATION', 'MOLWEIGHT', 'VOLUME', 'TID', 'TARGET_MAPPING', 'PUBLISHED_TYPE', 
        'PARENT_CMPD_CHEMBLID', 'CMPD_CHEMBLID', 'ACTIVITY_COMMENT', 'DOC_ID', 'PROTEIN_ACCESSION', 'BAO_FORMAT', 'CELL_CHEMBL_ID',
        'DATA_VALIDITY_COMMENT', 'CONFIDENCE_SCORE', 'ASSAY_SRC_DESCRIPTION', 'ASSAY_SRC_ID', 'BAO_ENDPOINT', 'PUBLISHED_UNITS', 
        'FIRST_PAGE', 'UO_UNITS', 'APD_CONFIDENCE', 'ALOGP', 'ASSAY_CHEMBLID', 'CANONICAL_SMILES', 'ASSAY_STRAIN', 'ASSAY_TAX_ID', 
        'MOL_PREF_NAME', 'PUBMED_ID', 'STANDARD_TYPE', 'CURATED_BY', 'STANDARD_VALUE', 'QUDT_UNITS', 'POTENTIAL_DUPLICATE', 'PARENT_MOLREGNO',
        'RELATION', 'STANDARD_UNITS', 'YEAR', 'DOC_CHEMBLID', 'ISSUE', 'ASSAY_ORGANISM', 'NUM_RO5_VIOLATIONS', 'APD_NAME', 'PUBLISHED_VALUE', 
        'DESCRIPTION', 'TARGET_CHEMBLID', 'JOURNAL', 'PSA', 'TARGET_TYPE', 'PREF_NAME', 'PCHEMBL_VALUE', 'COMPOUND_KEY', 'MOLREGNO', 'ASSAY_ID',
        'ORGANISM', 'ASSAY_TYPE']
        
        lcol_select = ["CMPD_CHEMBLID", "CANONICAL_SMILES", "STANDARD_VALUE", "PUBLISHED_UNITS", "STANDARD_TYPE", "PCHEMBL_VALUE", "CONFIDENCE_SCORE", "CURATED_BY", 
        "ASSAY_ID", "ASSAY_TYPE", "ASSAY_ORGANISM", "ASSAY_SRC_DESCRIPTION", "BAO_FORMAT", "BAO_ENDPOINT"]

        print "Size table: ", len(self.table)

        filout = open(pfilout, "w")
        filout.write("\t".join(lcol_select) + "\n")
        for row in self.table:
            lw = [row[h] for h in lcol_select]
            filout.write("\t".join(lw) + "\n")
        filout.close()

        self.pdataset = pfilout

    def writeTableAff(self, kaff = "PCHEMBL_VALUE"):

        pfilout = self.pr_out + "chem_aff_" + "-".join(self.ltypeAff)

        if path.exists(pfilout):
            self.paff = pfilout
            return pfilout

        filout = open(pfilout, "w")
        filout.write("CHEMBLID\tAff\tType\n")
        for row in self.table:
            filout.write(str(row["CMPD_CHEMBLID"]) + "\t" + str(row[kaff]) + "\t" + str(row["STANDARD_TYPE"]) + "\n")
        filout.close()

        self.paff = pfilout
        return pfilout



    def analysisTable(self, pranalysis):

        if not "paff" in self.__dict__:
            self.writeTableAff(pranalysis + "affAnalysis")

        if not "pdataset" in self.__dict__:
            self.CleanCHEMBLFileProtAff()

        # hist for activity
        runExternalSoft.histAffinity(self.paff)



    def checkIdenticSMI(self):

        lsmi = []
        i = 0
        imax = len(self.table)
        while i < imax:
            smi = self.table[i]["CANONICAL_SMILES"]
            smiclean = ligand.standardizeSMILES(smi)
            print "***", i, smiclean, "***"
            if smiclean == 1:
                i += 1
                continue

            if not smiclean in lsmi:
                lsmi.append(smiclean)
            else:
                del self.table[i]
                imax = imax - 1
                continue

            i += 1


    ##########################################
    # main for filtering with prot in function of Aff type and filin
    ##########################################

    def CleanCHEMBLFileProtAff(self):

        # if filout already exist
        pfilout = self.pr_out + "tab_filtered_" + "-".join(self.ltypeAff) + ".csv"

        if path.exists(pfilout):
            ltable = toolbox.matrixToList(pfilout)
            print len(ltable), "Nb selected compounds"
            print ltable[0].keys()
            self.table = ltable
            return ltable

        self.parseCHEMBLFile()
        print len(self.table), "Init cleaning"

        self.selectConfidencecore(cutoff=9)
        print len(self.table), "prot confidence"

        self.getOnlyExactConstant()
        print len(self.table), "strict value"

        self.getByTypeOfAff(self.ltypeAff)
        print len(self.table), self.ltypeAff

        self.MergeIdenticCHEMBLIDforACtivity()
        print len(self.table), "Repetition"

        self.selectAssayType("B")
        print len(self.table), "Type assay"

        # remove some biassay
        self.removeBA(self.outlier)
        print len(self.table), "remove bioassay"

        self.checkIdenticSMI()
        print len(self.table), "Identic SMI"

        self.writeTable(pfilout)

        # return the table in list
        return self.table



    def getChemSMI(self, pr_SMI_out):

        if not "table" in self.__dict__:
            print "ERROR: clean dataset first - run CleanCHEMBLFileProtAff"
        
        for chem in self.table:
            CHEMBLID = chem["CMPD_CHEMBLID"]
            pfilout = pr_SMI_out + CHEMBLID + ".smi"
            filout = open(pfilout, "w")
            filout.write(chem["CANONICAL_SMILES"])
            filout.close()




    def CleanCHEMBLFileCellLine(self):

        # add short cut if filtered table exist !!!!!

        self.parseCHEMBLFile()
        print len(self.table), "Init"

        self.getOnlyExactConstant()
        print len(self.table), "strict value"

        self.selectAssayType("F")
        print len(self.table), "Type assay"

        self.MergeIdenticCHEMBLIDforACtivity()
        print len(self.table), "Repetition"

        self.getByTypeOfAff(ltypeaff)
        print len(self.table), ltypeaff

        self.writeTable(pfilout)


        return self.table


