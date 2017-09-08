from os import listdir
from re import search
from shutil import copy2
from numpy import std, mean

import pathFolder
import runExternalSoft
import PDB
import parseShaep
import FPI
import parseSDF

class AnalyseClusterCpd:

    def __init__(self, sdata, pfilecluster, proutcluster, lprdockingpose):

        self.pcluster = pfilecluster
        self.prout = proutcluster
        self.lprdockingpose = lprdockingpose
        self.sdata = sdata

        # affinity
        if sdata != {}:
            daff = {}
            for compound in self.sdata:
                daff[compound["CMPD_CHEMBLID"]] = compound["PCHEMBL_VALUE"]
            self.Aff = daff

        prout = self.prout + "-".join(pfilecluster[0:-4].split("_")[1:]) + "/"
        prout = prout.replace(".", "")
        pathFolder.createFolder(prout)
        self.prout = prout

        filecluster = open(pfilecluster, "r")
        llinesCluster = filecluster.readlines()
        filecluster.close()
        dcluster = {}

        for lineCluster in llinesCluster[1:]:
            cluster = lineCluster.strip().split(",")[-1].replace("\"", "")
            #print cluster
            if not cluster in dcluster.keys():
                dcluster[cluster] = []
            compoundID = lineCluster.strip().split("\"")[1]
            dcluster[cluster].append(compoundID)

        self.clusters = dcluster


    def superimposedPoseCluster(self):

        for prdocking in self.lprdockingpose:
            prout = self.prout + prdocking.split("/")[-2] + "/"
            lposes = listdir(prdocking)
            proutSUperimpose = prout + "Superimposed/"
            pathFolder.createFolder(proutSUperimpose)
            for cluster in self.clusters.keys():
                pclusterpose = proutSUperimpose + str(cluster) + ".pdb"
                for compound in self.clusters[cluster]:
                    for pose in lposes:
                        if pose[-3:] != "sdf":
                            continue
                        elif search(compound, pose):
                            pposePDB = runExternalSoft.babelConvertSDFtoPDB(prdocking + pose)
                            cpose = PDB.PDB(PDB_input=pposePDB)
                            cpose.renameAtom()
                            cpose.writePDB(pfilout=pclusterpose, conect=1, model=1)



    def ShaepMatrix(self):

        d_sheap = {}
        for prdockingpose in self.lprdockingpose:
            lposesAll = listdir(prdockingpose)
            prout =  self.prout + prdockingpose.split("/")[-1]
            proutShaep = prout + "Sheap/"
            pathFolder.createFolder(proutShaep)
            for cluster in self.clusters.keys():
                lposecluster = []
                for compound in self.clusters[cluster]:
                    for poseAll in lposesAll:
                        pposetemp = prdockingpose + poseAll
                        if search(compound, poseAll) and search("pdb", poseAll):
                            if not pposetemp in lposecluster:
                                lposecluster.append(pposetemp)

                i = 0
                nb_pose = len(lposecluster)
                while i < nb_pose:
                    j = i
                    if not lposecluster[i] in d_sheap.keys():
                        d_sheap[lposecluster[i]] = {}
                    while j < nb_pose:
                        if not lposecluster[j] in d_sheap.keys():
                            d_sheap[lposecluster[j]] = {}
                        prtemp = prout + "tempSheap/"
                        pathFolder.createFolder(prtemp, clean=1)

                        # move file
                        pposei = prtemp + lposecluster[i].split("/")[-1]
                        pposej = prtemp + lposecluster[j].split("/")[-1]

                        copy2(lposecluster[i], pposei)
                        copy2(lposecluster[j], pposej)

                        psheap = runExternalSoft.runShaep(pposei, pposej, prtemp + "out", clean=1)
                        lsheapScore = parseShaep.parseOutputShaep(psheap)
                        d_sheap[lposecluster[i]][lposecluster[j]] = lsheapScore
                        d_sheap[lposecluster[j]][lposecluster[i]] = lsheapScore
                        j += 1
                    i += 1

                # write matrix
                for typescore in lsheapScore.keys():
                    pfiloutmatrix = proutShaep + str(typescore) + "_" + str(cluster)
                    filoutmatrix = open(pfiloutmatrix, "w")
                    lposes = [i.split("/")[-1].split(".")[0] for i in lposecluster]
                    filoutmatrix.write("\t".join(lposes)+"\n")
                    for posename in lposecluster:
                        filoutmatrix.write(posename.split("/")[-1].split(".")[0] + "\t")
                        lval = [str(d_sheap[posename][i][typescore]) for i in lposecluster]
                        filoutmatrix.write("\t".join(lval) + "\n")
                    filoutmatrix.close()
                    # affinity take in MDS folder !!!!
                    runExternalSoft.MatrixMCS(pfiloutmatrix, prout + "/MCS/" + str(cluster) + "_aff", pfiloutmatrix, 1)


    def FPIbycluster(self, pprot):


        d_FPI = {}

        prFPI = self.prout + "FPI/"
        pathFolder.createFolder(prFPI, clean=0)
        prFPItemp = self.prout + "FPItemp/"
        pathFolder.createFolder(prFPItemp, clean=1)


        lposesAll = listdir(self.prdocking)
        proutFPI = self.prout + "FPI/"
        pathFolder.createFolder(proutFPI)

        pfiloutaveragecluster = proutFPI + "meansInteraction"
        filoutaveragecluster = open(pfiloutaveragecluster, "w")
        lresglobal = []

        for cluster in self.clusters.keys():
            if not cluster in d_FPI.keys():
                d_FPI[cluster] = {}
            lposecluster = []
            for compound in self.clusters[cluster]:
                for poseAll in lposesAll:
                    pposetemp = self.prdocking + poseAll
                    if search(compound, poseAll) and search("pdb", poseAll):
                        if not pposetemp in lposecluster:
                            lposecluster.append(pposetemp)

            # add pose in protein file
            lResBS = []
            for posecluster in lposecluster:
                cProt = PDB.PDB(pprot, hydrogen=1)
                cProt.addLigand(posecluster)
                CFPI = FPI.ligFPI(cProt, prFPItemp, ligID="LIG")
                CFPI.computeFPI(clean=1)
                print CFPI
                d_FPI[cluster][posecluster] = CFPI

                ligID = CFPI.FPI.keys()[0]
                for resBS in CFPI.FPI[ligID]:
                    if not resBS in lResBS:
                        lResBS.append(resBS)
                    if not resBS in lresglobal:
                        lresglobal.append(resBS)

            # write out put
            filout = open(proutFPI + "fpi_out_" + str(cluster), "w")
            filout.write("pose\t" + " ".join(lResBS) + "\n")
            for pose in lposecluster:
                filout.write(pose.split("/")[-1])
                for resBS in lResBS:
                    if not resBS in d_FPI[cluster][pose].FPI[ligID].keys():
                        filout.write("\t0000000")
                    else:
                        filout.write("\t" + str(d_FPI[cluster][pose].FPI[ligID][resBS]))
                filout.write("\n")
            filout.close()


        # average interaction
        filoutaveragecluster.write("Clusters\t" + "\t".join(lresglobal) + "\n")
        for cluster in d_FPI.keys():
            filoutaveragecluster.write(str(cluster))

            for resglobal in lresglobal:
                lbit = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                #
                for pose in d_FPI[cluster].keys():
                    if resglobal in d_FPI[cluster][pose].FPI[ligID].keys():
                        for i in range(0, 7):
                            print lbit[i] , d_FPI[cluster][pose].FPI[ligID][resglobal][i]
                            lbit[i] = float(lbit[i]) + float(d_FPI[cluster][pose].FPI[ligID][resglobal][i])
                        for i in range(0, 7):
                            lbit[i] = str(lbit[i] / len(d_FPI[cluster].keys()))
                    else:
                        lbit = [str(i) for i in lbit]
                filoutaveragecluster.write("\t" + " ".join(lbit))
            filoutaveragecluster.write("\n")
        filoutaveragecluster.close()


    def summarize(self, opaff = 1):

        for prdockingpose in self.lprdockingpose:
            prout = self.prout + prdockingpose.split("/")[-2] + "/"
            prsum = prout + "Summary/"
            print prsum
            pathFolder.createFolder(prsum)
            lposesAll = listdir(prdockingpose)

            pfilout = prsum + "summary"
            filout = open(pfilout, "w")
            if opaff == 0:
                filout.write("ID\tMDockingScore\tSDDockingScore\tMEmodel\tSDEmodel\tNB\n")
            else:
                filout.write("ID\tMDockingScore\tSDDockingScore\tMEmodel\tSDEmodel\tMAff\tSDAff\tNb\n")
            for cluster in self.clusters.keys():
                ldockingScore = []
                lemodel = []
                laff = []
                for compound in self.clusters[cluster]:
                    for poseAll in lposesAll:
                        if search(compound + ".1.sdf", poseAll):
                            pposetemp = prdockingpose + poseAll
                            break
                    if pposetemp == "": continue
                    csdf = parseSDF.sdf(pposetemp)
                    ddock = csdf.get_dockingscore()
                    pposetemp = ""
                    ldockingScore.append(ddock[compound]["r_i_docking_score"])
                    lemodel.append(ddock[compound]["r_i_glide_emodel"])
                    if opaff != 0: laff.append(float(self.Aff[compound]))

                MdockingScore = mean(ldockingScore)
                SDdockingscore = std(ldockingScore)

                Memodel = mean(lemodel)
                SDemodel = std(lemodel)
                nb = len(lemodel)

                if opaff == 1:
                    Maff = mean(laff)
                    SDaff = std(laff)
                    filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (cluster, MdockingScore, SDdockingscore, Memodel, SDemodel, Maff, SDaff, nb))
                else:
                    filout.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (cluster, MdockingScore, SDdockingscore, Memodel, SDemodel, nb))
            filout.close()

