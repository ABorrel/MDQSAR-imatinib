from copy import deepcopy
from numpy import mean, std
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
            else:
                i += 1

        print len(self.table)

    def getOnlyIC50(self):

        if not "table" in dir(self):
            self.parseCHEMBLFile()

        i = 0
        imax = len(self.table)
        while i < imax:
            row = self.table[i]
            if row["PUBLISHED_TYPE"] != "IC50":
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
            # case not problem
            #print d_CHEMBLID.keys()[i]
            if len(d_CHEMBLID[d_CHEMBLID.keys()[i]]) == 1:
                i += 1
                continue
            else:
                # Control the published value
                l_PUBLISHED_VALUE = [float(d_CHEMBLID[d_CHEMBLID.keys()[i]][k]["PUBLISHED_VALUE"]) for k in range(0, len(d_CHEMBLID[d_CHEMBLID.keys()[i]]))]
                l_PUBLISHED_UNITS = [d_CHEMBLID[d_CHEMBLID.keys()[i]][k]["PUBLISHED_UNITS"] for k in range(0, len(d_CHEMBLID[d_CHEMBLID.keys()[i]]))]

                if not l_PUBLISHED_UNITS.count(l_PUBLISHED_UNITS[0]) == len(l_PUBLISHED_UNITS):
                    l_PUBLISHED_VALUE = toolbox.convertUnit(l_PUBLISHED_VALUE, l_PUBLISHED_UNITS)
                    #print l_PUBLISHED_UNITS

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




