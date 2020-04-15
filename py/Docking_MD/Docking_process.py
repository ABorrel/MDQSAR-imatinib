from PIL import Image, ImageFont, ImageDraw
import pathFolder
from os import path
import toolbox
import runExternalSoft
#from shutil import copy

#import runExternalSoft
#import pathFolder
#import RMSD

# import parsor scripts
import sys
sys.path.insert(0, "./Parser/") # for window dev
import parseSDF

font = ImageFont.truetype("./OpenSans-Regular.ttf", size=24)

class Docking_process:
    def __init__(self, psdf_poses, pr_out):

        self.pr_out = pr_out
        self.psdf_poses = psdf_poses

    def loadSDF(self, splitSDF=0):
        """
        Load the pose results
        """
        self.cPoses = parseSDF.parseSDF(self.pr_out)
        self.cPoses.parseSDFFile(self.psdf_poses)
        
        if splitSDF == 1:
            pr_poses = pathFolder.createFolder(self.pr_out + "ALL_POSES/")
            self.cPoses.splitChem(pr_poses)


    def get_bestdockingscore(self):
        """Keep smaller score by chem"""

        if not "lc" in self.cPoses.__dict__:
            print "Load first sdf with poses"
            return 

        if "dscores" in self.__dict__:
            return self.dscores

        dscores = {}
        for dchem in self.cPoses.lc:
            # case where protein is included, case of XP docking
            if not "r_i_docking_score" in dchem.keys():
                continue

            chemblID = dchem["s_m_entry_name"].split(".")[0]
            #print chemblID

            if not chemblID in dscores.keys():
                dscores[chemblID] = {}
                dscores[chemblID]["count"] = 1
            else:
                dscores[chemblID]["count"] = dscores[chemblID]["count"] + 1

            if not "r_i_docking_score" in dscores[chemblID].keys():
                dscores[chemblID]["r_i_docking_score"] = float(dchem["r_i_docking_score"])
                dscores[chemblID]["r_i_glide_emodel"] = float(dchem["r_i_glide_emodel"])
            else:
                if float(dchem["r_i_docking_score"]) < dscores[chemblID]["r_i_docking_score"]:
                    dscores[chemblID]["r_i_docking_score"] = float(chemblID["r_i_docking_score"])
                    dscores[chemblID]["r_i_glide_emodel"] = float(chemblID["r_i_glide_emodel"])

        self.dscores = dscores

        # write
        pfilout = self.pr_out + "score_poses.txt"
        filout = open(pfilout, "w")
        filout.write("Chemicals\tNb poses\tGlide score\temodel score\n")
        for chemblID in dscores.keys():
            filout.write("%s\t%s\t%s\t%s\n"%(chemblID, dscores[chemblID]["count"], dscores[chemblID]["r_i_docking_score"], dscores[chemblID]["r_i_glide_emodel"]))
        filout.close()


    def plot_dockScoreVSActivity(self, ptableCHEMBL):

        pfilout = self.pr_out + "ScoreVSAff.txt"
        if path.exists(pfilout):
            return pfilout
        
        # load ChEMBL dataset
        ldchem = toolbox.matrixToList(ptableCHEMBL)


        filout = open(pfilout, "w")
        filout.write("IDCHEMBL\tDock_score\temodel_score\tAff\ttypeAff\tNB poses\n")

        self.dscores

        for dchem in ldchem:
            print(chemID)
            filout.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(dchem["CMPD_CHEMBLID"], self.dscores[chemID]["r_i_docking_score"], self.dscores[chemID]["r_i_glide_emodel"], dchem[chemID]["PCHEMBL_VALUE"], dchem[chemID]["STANDARD_TYPE"], self.dscores[chemID]["count"]))
        
        filout.close()
        runExternalSoft.corPlot(pfilout, ptableCHEMBL, prout)



    def get_bestPoses(self):
        """
        Function used to extract best poses based on the min docking score
        """
        
        if not "dscores" in self.__dict__:
            self.get_bestdockingscore()
        
        # define directory
        pr_bestPoses = pathFolder.createFolder(self.pr_out + "BEST_POSES/")
        self.pr_bestPoses = pr_bestPoses
        lbest_poses = []

        for chemID in self.dscores.keys():
            bestScoreDock = float(self.dscores[chemID]["r_i_docking_score"])
            bestScoreEmodel = float(self.dscores[chemID]["r_i_glide_emodel"])

            ipose = 0
            nbpose = len(self.cPoses.lc)
            while ipose < nbpose:
                namepose = self.cPoses.lc[ipose]["s_m_entry_name"]
                nameChem = namepose.split(".")[0]
                try:
                    dock_score = float(self.cPoses.lc[ipose]["r_i_docking_score"])
                    emodel_score = float(self.cPoses.lc[ipose]["r_i_glide_emodel"])
                except:
                    ipose = ipose + 1
                    continue
                if chemID == nameChem:
                    if bestScoreDock == dock_score and bestScoreEmodel == emodel_score:
                        pfilout = self.pr_bestPoses + nameChem + ".sdf"
                        if not path.exists(pfilout) or path.getsize(pfilout) < 20: # control if exist
                            filout = open(pfilout, "w")
                            filout.write(self.cPoses.lc[ipose]["sdf"])
                            filout.close()
                            # apply a format with babel to have a proper sdf
                            runExternalSoft.babelConverttoSDF(pfilout)
                            break

                    else:
                        lbest_poses.append(self.cPoses.lc[ipose])
                ipose = ipose + 1

        print "NB best poses:", len(self.cPoses.lc[ipose])


    def get_TopRanking(self, ptableCHEMBL, nrank = 5):

        pr_topChem = pathFolder.createFolder(self.pr_out + "top-" + str(nrank) + "/")
        
        # load dock score
        if not "dscores" in self.__dict__:
            self.get_bestdockingscore()

        # load dataset
        ldchem = toolbox.matrixToList(ptableCHEMBL)
        
        dbytypeAct = {}
        for dchem in ldchem:
            typeAff = dchem["STANDARD_TYPE"]
            if not typeAff in dbytypeAct.keys():
                dbytypeAct[typeAff] = {}
            nameChem = dchem["CMPD_CHEMBLID"]
            dbytypeAct[typeAff][nameChem] = dchem

        # extract top 10 dock score
        for typeAff in dbytypeAct.keys():
            print typeAff
            laff = []
            for pose in self.cPoses.lc:
                pose_chem = pose["s_m_entry_name"].split(".")[0]
                if pose_chem in dbytypeAct[typeAff].keys():
                    laff.append(float(pose["r_i_docking_score"]))

            # order
            laff.sort()
            i = 0
            irank = 1
            l_chem_image = []
            while irank <= nrank:
                for pose in self.cPoses.lc:
                    pose_chem = pose["s_m_entry_name"].split(".")[0]
                    try:score_dock = float(pose["r_i_docking_score"])
                    except:continue
                    if score_dock == laff[i]:
                        if pose_chem in l_chem_image:
                            continue
                        else:
                            l_chem_image.append(pose_chem)
                        # path png and sdf
                        p_image = "%s%s_%s-%s.png"%(pr_topChem, irank, pose_chem, typeAff) 
                        p_smiout = "%s%s_%s-%s.smi"%(pr_topChem, irank, pose_chem, typeAff) 

                        # to write on figure
                        lw = ["ChEMBL ID: " + str(pose_chem)]
                        lw.append("Rank: " + str(irank))
                        lw.append("p" + typeAff + ": " + str(dbytypeAct[typeAff][pose_chem]["PCHEMBL_VALUE"]))
                        lw.append("Docking score: {:0.2f}".format(float(pose["r_i_docking_score"])))
                        lw.append("Emodel score: {:0.2f}".format(float(pose["r_i_glide_emodel"])))

                        # generate png
                        f_smiout = open(p_smiout, "w")
                        f_smiout.write(dbytypeAct[typeAff][pose_chem]["CANONICAL_SMILES"])
                        f_smiout.close()

                        # png
                        runExternalSoft.molconvert(p_smiout, p_image)

                        img = Image.open(p_image)
                        imgnew = Image.new("RGBA", (580, 775), (250, 250, 250))
                        imgnew.paste(img, (0,0))
                        draw = ImageDraw.Draw(imgnew)
                        draw.text((10, 600), lw[0], (0, 0, 0), font=font)
                        draw.text((10, 625), lw[1], (0, 0, 0), font=font)
                        draw.text((10, 650), lw[2], (0, 0, 0), font=font)
                        draw.text((10, 675), lw[3], (0, 0, 0), font=font)
                        draw.text((10, 700), lw[4], (0, 0, 0), font=font)
                        imgnew.save(p_image)
                        irank = irank + 1
                        break
                
                i = i + 1









"""
def plotRMSDVSDockingScore(ddockingscore, ltabCHEMBL, ptableCHEMBL,prMDanalysis, prout):


    pfilout = prout + "ScoreVSRMSD"
    filout = open(pfilout, "w")
    filout.write("IDCHEMBL\tDock_score\temodel\tRMSDca\tRMSDall\tDmax\tRMSDlig\tAff\ttypeAff\n")

    for daff in ltabCHEMBL:
        CHEMBLid = daff["CMPD_CHEMBLID"]
        pRMSDin = prMDanalysis + CHEMBLid + "_2hyy_MD/RMSDs/"
        if not path.exists(pRMSDin):
            continue

        RMSDChem = RMSD.RMSD(pRMSDin)
        RMSDChem.loadRMSDs(["ligand", "protein"])
        RMSDprot = RMSDChem.MRMSDprot()
        RMSDlig = RMSDChem.MRMSDlig()

        print RMSDprot, RMSDlig


        filout.write("%s\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\n"%(CHEMBLid,
                                                                 ddockingscore[CHEMBLid]["r_i_docking_score"],
                                                                 ddockingscore[CHEMBLid]["r_i_glide_emodel"],
                                                                 RMSDprot[0], RMSDprot[1], RMSDprot[2], RMSDlig,
                                                                 daff["PCHEMBL_VALUE"], daff["STANDARD_TYPE"]))


    filout.close()
    runExternalSoft.corPlot(pfilout, ptableCHEMBL, prout, "RMSD")

"""