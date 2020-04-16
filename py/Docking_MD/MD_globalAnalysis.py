from os import path, listdir
from re import search
import toolbox
import calculate
import runExternalSoft
import pathFolder

import sys
sys.path.insert(0, "./Parser/")
import PDB
import parseTMalign



# class to analyse MD after computed MD and RMSD
class MD_globalAnalysis:
    def __init__(self, p_dataset, pr_MDrun, pr_MDout, pr_out):
        self.p_dataset = p_dataset
        self.pr_MDout = pr_MDout
        self.pr_MDrun = pr_MDrun
        self.pr_out = pr_out

    def buildRMDSSheets(self):
        """
        Function use to build a pdf by chemicals with RMSD prot, ligand, BS and RMSF by residues of protein
        Need 3 files, RMSD for prot, RMSD for ligand and RMSF
        """
        l_fMDout = listdir(self.pr_MDout)

        # load ChEMBL table
        ldchem = toolbox.matrixToList(self.p_dataset)
        for dchem in ldchem:
            dfile = {}
            dfile["prot RMSD"] = ""
            dfile["lig RMSD"] = ""
            dfile["RMSF residue"] = ""

            ChEMBL_id = dchem["CMPD_CHEMBLID"]
            typeAff = dchem["STANDARD_TYPE"]
            # extract RMSD files
            for fMDout in l_fMDout:
                ChEMBL_id_folder = fMDout.split("_")[0]
                if ChEMBL_id == ChEMBL_id_folder:
                    
                    
                    # protein
                    p_RMSD_prot = self.pr_MDout + fMDout + "/RMSDs/protein/protRMSD"
                    if path.exists(p_RMSD_prot):
                        dfile["prot RMSD"] = p_RMSD_prot
                    else:
                        print "Error in " + fMDout + ": prot RMSD missing"
                        self.computeRMSDProt(self.pr_MDout + fMDout + "/")

                    # lig need to define here
                    p_RMSD_lig = self.pr_MDout + fMDout + "/RMSDs/ligand/ligRMSD"
                    if path.exists(p_RMSD_lig):
                        dfile["lig RMSD"] = p_RMSD_lig
                    else:
                        print "Error in " + fMDout + ": lig RMSD missing"
                        p_RMSD_lig = self.computeRMSDLig(self.pr_MDout + fMDout + "/")
                        dfile["lig RMSD"] = p_RMSD_lig

                    
                    # RMSF with the binding site represented
                    p_RMSF_res = self.pr_MDout + fMDout + "/RMSDs/residues/resRMSD_BS"
                    if path.exists(p_RMSF_res):
                        dfile["RMSF residue"] = p_RMSF_res
                    else:
                        print "Error in " + fMDout + ": lig RMSF BS missing"
                        p_RMSF_res = self.computeRMSFresBS(self.pr_MDout + fMDout + "/")
                        dfile["lig RMSD"] = p_RMSF_res
                
                    # build figure 3 panels
                    runExternalSoft.RMSD3panels(dfile["prot RMSD"], dfile["lig RMSD"], dfile["RMSF residue"], ChEMBL_id, pathFolder.createFolder(self.pr_out + typeAff + "/"))


    def computeRMSDProt(self, pr_MDout):
         # load ligand in frame 0
        cfram0 = PDB.PDB(pr_MDout + "framesMD/frame_00000.pdb")
        cfram0.get_lAtoms()

        pr_TMalign = pr_MDout + "RMSDs/superimpose/"
        l_pTMaling = listdir(pr_TMalign)

        dRMSD = {}
        for pTMalign in l_pTMaling:
            frame = pTMalign.split("_")[-1]
            #print frame
            dmatrixTMalign = toolbox.loadMatrixTMalign(pr_TMalign + pTMalign)
            cFrame = PDB.PDB( "%sframesMD/frame_%s.pdb"%(pr_MDout, frame))
            cFrame.get_lAtoms()

            for atomLig in cFrame.latom:
                atomLig.applyMatrixRotTransloc(dmatrixTMalign)
            
            RMSDframe = calculate.RMSDTwoList(cfram0.latom, cFrame.latom)
            dRMSD[frame] = RMSDframe
        
        # write the RMSD lig file
        pfilout = pr_MDout + "RMSDs/protein/protRMSD_all"
        filout = open(pfilout, "w")
        filout.write("Time\tRMSDall\tRMSDC\tDmax\n")
        filout.write("0.0\t0.0\t0.0\t0.0\n")

        i = 1
        imax = len(dRMSD.keys())
        while i <= imax:
            frame = str("%05d" % (i))
            filout.write("%.2f\t%s\t%s\t%s\n"%(i/100.0, dRMSD[frame][0], dRMSD[frame][1], dRMSD[frame][2]))
            i = i + 1
        filout.close()
        return pfilout



    def computeRMSDLig(self, pr_MDout):
        
        # load ligand in frame 0
        clig0 = PDB.PDB(pr_MDout + "lig/LGD_00000.pdb")
        clig0.get_lAtoms()

        pr_TMalign = pr_MDout + "RMSDs/superimpose/"
        l_pTMaling = listdir(pr_TMalign)

        dRMSD = {}
        for pTMalign in l_pTMaling:
            frame = pTMalign.split("_")[-1]
            #print frame
            dmatrixTMalign = toolbox.loadMatrixTMalign(pr_TMalign + pTMalign)
            cligFrame = PDB.PDB( "%slig/LGD_%s.pdb"%(pr_MDout, frame))
            cligFrame.get_lAtoms()

            for atomLig in cligFrame.latom:
                atomLig.applyMatrixRotTransloc(dmatrixTMalign)
            
            RMSDframe = calculate.RMSDTwoList(clig0.latom, cligFrame.latom)
            dRMSD[frame] = RMSDframe[0]
        
        # write the RMSD lig file
        pfilout = pr_MDout + "RMSDs/ligand/ligRMSD"
        filout = open(pfilout, "w")
        filout.write("Time\tRMSD\n")
        filout.write("0.0\t0.0\n")

        i = 1
        imax = len(dRMSD.keys())
        while i <= imax:
            frame = str("%05d" % (i))
            filout.write("%.2f\t%s\n"%(i/100.0, dRMSD[frame]))
            i = i + 1
        filout.close()
        return pfilout


    def computeRMSFresBS(self, pr_MDout):
        

        # load BS in frame 0
        l_pBS = listdir(pr_MDout + "BSs/")
        l_res = []
        for pBS in l_pBS:
            cBS = PDB.PDB(pr_MDout + "BSs/" + pBS)
            dres = cBS.get_byres()
            for res in dres.keys():
                nRes = res.split("_")[1]
                if not nRes in l_res:
                    l_res.append(nRes)

        # rewrite RMSF with binding site
        presRMSF = pr_MDout + "RMSDs/residues/resRMSD"
        ldresRMSF = toolbox.matrixToList(presRMSF)

        # rewrting 
        pfilout = pr_MDout + "RMSDs/residues/resRMSD_BS"
        filout = open(pfilout, "w")
        filout.write("NameRes\tall\tCa\tDmax\tBS\n")
        for dresRMSF in ldresRMSF:
            if dresRMSF["NameRes"] in l_res:
                BS = 1
            else:
                BS = 0
            
            filout.write("%s\t%s\t%s\t%s\t%s\n"%(dresRMSF["NameRes"], dresRMSF["all"], dresRMSF["Ca"], dresRMSF["Dmax"], BS))

        filout.close()
        return pfilout



    def mergeRMSDSheets(self):

        lfolder = listdir(self.pr_out)
        for folderAffType in lfolder:
            lpdfs = listdir(self.pr_out + folderAffType + "/")
            i = 0
            imax = len(lpdfs)
            while i < imax:
                lpdfs[i] = self.pr_out + folderAffType + "/" + lpdfs[i]
                i + i + 1
        
        pout = self.pr_out + "RMSD_" + folderAffType + ".pdf"
        runExternalSoft.mergepdfs(lpdfs, pout)

        return

