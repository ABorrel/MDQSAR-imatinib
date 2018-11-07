from PIL import Image, ImageFont, ImageDraw
from os import path
from shutil import copy

import runExternalSoft
import pathFolder
import RMSD

font = ImageFont.truetype("OpenSans-Regular.ttf", size=24)


def dockingScoreAnalysis(ddockingscore, ltabCHEMBL, ptableCHEMBL, prout):

    pfilout = prout + "ScoreVSAff.txt"
    if not path.exists(pfilout):
        filout = open(pfilout, "w")
        filout.write("IDCHEMBL\tDock_score\temodel\tAff\ttypeAff\n")

        for daff in ltabCHEMBL:
            try: filout.write(str(daff["CMPD_CHEMBLID"]) + "\t" + str(ddockingscore[daff["CMPD_CHEMBLID"]]["r_i_docking_score"])
                              + "\t" + str(ddockingscore[daff["CMPD_CHEMBLID"]]["r_i_glide_emodel"])
                              + "\t" + str(daff["PCHEMBL_VALUE"]) + "\t" + str(daff["STANDARD_TYPE"]) + "\n")
            except: pass
        filout.close()
    runExternalSoft.corPlot(pfilout, ptableCHEMBL, prout)



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




def rankingTop(spose, ltabCHEMBL, pbestposes, prout, nrank = 50):

    #order value
    daff = {}
    for daffCHEM in ltabCHEMBL:
        if not daffCHEM["STANDARD_TYPE"] in daff.keys():
            daff[daffCHEM["STANDARD_TYPE"]] = {}
        CHEMBLID = daffCHEM["CMPD_CHEMBLID"]
        daff[daffCHEM["STANDARD_TYPE"]][CHEMBLID] = float(daffCHEM["PCHEMBL_VALUE"])


    #for afftype in daff.keys():
    for afftype in ["Ki"]:
        lval = daff[afftype].values()
        lval.sort(reverse=True)

        lflag = []
        rank = 1
        for val in lval:
            for CHEMID in daff[afftype].keys():
                if daff[afftype][CHEMID] == val:
                    if not CHEMID in spose.docking.keys():
                        continue
                    if CHEMID in lflag:
                        continue
                    else:
                        lflag.append(CHEMID)
                    lw = ["ChEMBL ID: " + str(CHEMID)]
                    lw.append("Rank: " + str(rank))
                    lw.append("p" + afftype + ": " + str(val))
                    lw.append("Docking score: {:0.2f}".format(float(spose.docking[CHEMID]["r_i_docking_score"])))
                    lw.append("Emodel score: {:0.2f}".format(float(spose.docking[CHEMID]["r_i_glide_emodel"])))

                    pimageout = prout + str(rank) + "_" + str(CHEMID) + "_" + str(afftype) + ".png"
                    psdfout = prout + str(rank) + "_" + str(CHEMID) + "_" + str(afftype) + ".sdf"
                    pCHEMBLpng = pathFolder.PR_PNG + CHEMID + ".png"
                    psdf = pbestposes + CHEMID + ".1.sdf"

                    if not path.exists(pCHEMBLpng) or path.getsize(pCHEMBLpng) < 10:
                        print "not found", pCHEMBLpng
                        rank = rank + 1
                        continue
                    else:
                        img = Image.open(pCHEMBLpng)
                        imgnew = Image.new("RGBA", (580, 775), (250, 250, 250))
                        imgnew.paste(img, (0,0))
                        draw = ImageDraw.Draw(imgnew)
                        draw.text((10, 600), lw[0], (0, 0, 0), font=font)
                        draw.text((10, 625), lw[1], (0, 0, 0), font=font)
                        draw.text((10, 650), lw[2], (0, 0, 0), font=font)
                        draw.text((10, 675), lw[3], (0, 0, 0), font=font)
                        draw.text((10, 700), lw[4], (0, 0, 0), font=font)
                        imgnew.save(pimageout)

                        copy(psdf, psdfout)

                        rank = rank + 1

            if rank > nrank:
                break

