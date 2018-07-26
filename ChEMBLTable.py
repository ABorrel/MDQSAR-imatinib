from copy import deepcopy
from numpy import mean, std
from os import path

import collections
import runExternalSoft
import toolbox

class CHEMBL:
    def __init__(self, pfilin):
        self.pfilin = pfilin

    def parseCHEMBLFile(self):

        lout = []
        filin = open(self.pfilin, "r")
        llines = filin.readlines()
        filin.close()

        lhead = llines[0].split("\t")[0:-1]
        print lhead

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

        print len(self.table)

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
            print i, d_CHEMBLID.keys()[i]
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

        print len(self.table)
        print self.table[0]

        filout = open(pfilout, "w")
        lheader = self.table[0].keys()
        filout.write("\t".join(lheader) + "\n")
        for row in self.table:
            lw = [row[h] for h in lheader]
            filout.write("\t".join(lw) + "\n")


        filout.close()


    def writeTableAff(self, pfilout, kaff = "PCHEMBL_VALUE"):

        if path.exists(pfilout):
            return pfilout

        filout = open(pfilout, "w")
        filout.write("CHEMBLID\tAff\tType\n")
        for row in self.table:
            filout.write(str(row["CMPD_CHEMBLID"]) + "\t" + str(row[kaff]) + "\t" + str(row["STANDARD_TYPE"]) + "\n")
        filout.close()

        self.paff = pfilout
        return pfilout



    def analysisTable(self, pranalysis):

        if not "paff" in dir(self):
            self.writeTableAff(pranalysis + "affAnalysis")


        runExternalSoft.histAffinity(self.paff)
