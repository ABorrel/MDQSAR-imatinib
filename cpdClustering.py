from os import listdir
from re import search

import pathFolder
import runExternalSoft
import PDB

class AnalyseClusterCpd:

    def __init__(self, pfilecluster, proutcluster, prdockingpose):

        self.pcluster = pfilecluster
        self.prout = proutcluster
        self.prdocking = prdockingpose

        prout = self.prout + "-".join(pfilecluster[0:-4].split("_")[1:]) + "/"
        pathFolder.createFolder(prout)
        self.prout = prout

        filecluster = open(pfilecluster, "r")
        llinesCluster = filecluster.readlines()
        filecluster.close()
        dcluster = {}

        for lineCluster in llinesCluster[1:]:
            cluster = lineCluster.strip().split(",")[-1].replace("\"", "")
            print cluster
            if not cluster in dcluster.keys():
                dcluster[cluster] = []
            compoundID = lineCluster.strip().split("\"")[1]
            dcluster[cluster].append(compoundID)

        self.clusters = dcluster


    def superimposedPoseCluster(self):

        lposes = listdir(self.prdocking)
        for cluster in self.clusters.keys():
            pclusterpose = self.prout + str(cluster) + ".pdb"
            for compound in self.clusters[cluster]:
                for pose in lposes:
                    if search(compound, pose):
                        pposePDB = runExternalSoft.babelConvertSDFtoPDB(self.prdocking + pose)
                        cpose = PDB.PDB(PDB_input=pposePDB)
                        cpose.writePDB(pfilout=pclusterpose, model=1)



