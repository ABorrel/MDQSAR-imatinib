#!/usr/bin/env python
# -*- coding:utf-8 -*-
from time import time
import getopt, sys, os, copy, glob, re
from openbabel import OBMol, OBConversion
from optparse import OptionParser
from bitarray import bitarray

from ring import *
from interactions import *
import runExternalSoft

def get_FPI(pligPDB, ppocketPDB, lres, filout):

    x = time()

    # convert PDB to mol2
    pligmol2 = runExternalSoft.babelPDBtoMOL2(pligPDB)
    ppocketMol2 = runExternalSoft.babelPDBtoMOL2(ppocketPDB)

    # opening the molecule files
    conv = OBConversion()

    # remove because format change
    conv.SetInFormat("mol2")

    st_protref = OBMol()
    st_ligref = OBMol()

    conv.ReadFile(st_protref, ppocketMol2)
    conv.ReadFile(st_ligref, pligmol2)

    # creating backbone atom list
    protrefhandle = open(ppocketMol2, 'r')
    protreflines = [line for line in protrefhandle]
    protrefhandle.close()

    l_atom_backbone = ['C', 'N', 'O', 'H', 'CA', 'HA']
    l_serial_backbone_ref = []

    for line in protreflines:
        if len(line.split()) > 1:
            if line.split()[1] in l_atom_backbone:
                l_serial_backbone_ref.append(int(line.split()[0]))

    # fingerprint reference
    d_ref_choice = getresiduedict(st_protref, lres)
    d_refring = getringdict(st_protref)
    ringinteraction(d_ref_choice, d_refring, lres, st_protref, st_ligref)
    otherinteractions(d_ref_choice, lres, st_protref, st_ligref, l_serial_backbone_ref)
    hbonddockprot(d_ref_choice, lres, st_protref, st_ligref)

    # write reference
    filout.write(pligPDB.split("/")[-1][0:-4] + "\t" + "-".join(lres) + "\t")
    lbit = [str(collectbit(d_ref_choice, [res])) for res in lres]
    filout.write("-".join(lbit))
    filout.write("\n")

    y = time()
    print 'Total time taken %.3f s.' % (y - x)

    dout = {}
    for res in lres:
        dout[res] = str(collectbit(d_ref_choice, [res]))

    return dout


def collectbit(dict1, residuechoice):
    """This function will return bit array in string form"""

    bit = bitarray()
    for residue in residuechoice:
        bit.extend(dict1[residue])

    stringbit = str(bit).split("'")[1]

    return stringbit


"""
def FingerPrintInteraction(p_pocket_ref, p_lig_ref, pr_pose_and_pocket, l_res, p_filout, debug=1):
    ""
	update from -> Radifar, M., Yuniarti, N., & Istyastono, E. P. (2013). PyPLIF: Python-based Protein-Ligand Interaction Fingerprinting. Bioinformation, 9(6), 325–8. doi:10.6026/97320630009325
	change loader and remove config file
	But need fixe corretly charge => the ref using PLANTS 
	And why two different BS fix or not only one
	""

    #load file result if exist
    if os.path.exists(p_filout):
        if os.path.getsize(p_filout) > 400:
            return loader.LoadFPIResults(p_filout)

    x = time()
    filout = open(p_filout, "w")
    d_out = {}

    print "========"
    print p_pocket_ref
    print p_lig_ref
    print l_res
    print pr_pose_and_pocket
    print p_filout
    print "========"

    # define
    d_file = {}
    l_file_pose = os.listdir(pr_pose_and_pocket)
    for file_pose in l_file_pose:
        name_lig_pose = file_pose[0:-5].split("_")[0] + "_" + file_pose[0:-5].split("_")[1]
        if re.search("pocket.mol2", file_pose):
            if not name_lig_pose in d_file:
                d_file[name_lig_pose] = {}
            d_file[name_lig_pose]["pocket"] = pr_pose_and_pocket + file_pose

        elif re.search(".mol2", file_pose) and not re.search("protein", file_pose):
            if not name_lig_pose in d_file:
                d_file[name_lig_pose] = {}
            d_file[name_lig_pose]["pose"] = pr_pose_and_pocket + file_pose

    # opening the molecule files
    conv = OBConversion()

    # remove because format change
    conv.SetInFormat("mol2")

    st_protref = OBMol()
    st_ligref = OBMol()
    st_ligdock = OBMol()
    st_protdock = OBMol()

    conv.ReadFile(st_protref, p_pocket_ref)
    conv.ReadFile(st_ligref, p_lig_ref)

    # creating backbone atom list
    protrefhandle = open(p_pocket_ref, 'r')
    protreflines = [line for line in protrefhandle]
    protrefhandle.close()

    pbfhandle = open(p_pocket_ref, 'r')
    pbflines = [line for line in pbfhandle]
    pbfhandle.close()
    l_atom_backbone = ['C', 'N', 'O', 'H', 'CA', 'HA']
    l_serial_backbone_ref = []
    l_serial_backbone_test = []

    for line in protreflines:
        if len(line.split()) > 1:
            if line.split()[1] in l_atom_backbone:
                l_serial_backbone_ref.append(int(line.split()[0]))

    for line in pbflines:
        if len(line.split()) > 1:
            if line.split()[1] in l_atom_backbone:
                l_serial_backbone_test.append(int(line.split()[0]))

    # fingerprint reference
    d_ref_choice = getresiduedict(st_protref, l_res)
    d_refring = getringdict(st_protref)
    ringinteraction(d_ref_choice, d_refring, l_res, st_protref, st_ligref)
    otherinteractions(d_ref_choice, l_res, st_protref, st_ligref, l_serial_backbone_ref)
    hbonddockprot(d_ref_choice, l_res, st_protref, st_ligref)

    # bit reference
    # d_bit_ref = collectbit(d_ref_choice, l_res)

    # write reference
    filout.write("Name pose\tScore\t" + "\t".join(l_res) + "\n")
    filout.write("REF_" + p_lig_ref.split("/")[-1][0:-5] + "\t" + "1.0")
    for res in l_res:
        filout.write("\t" + str(collectbit(d_ref_choice, [res])))
    filout.write("\n")

    # implement output
    # 	d_out["FPI"] = collectbit(d_ref_choice, l_res)
    # 	d_out["score"] = 1.0

    # initialisation
    d_template = getresiduedict(st_protref, l_res)  # define template

    for name_lig_pose in d_file.keys():

        name_lig = name_lig_pose.split("_")[0]
        name_pose = name_lig_pose.split("_")[1]

        d_check_choice = copy.deepcopy(d_template)

        # print pose, d_file[name_lig_pose]
        p_pose = d_file[name_lig_pose]["pose"]
        p_pocket = d_file[name_lig_pose]["pocket"]

        if debug == 1:
            print "===================="
            print p_pose, "Path pose"
            print p_pocket, "Path pocket"
            print d_check_choice, "Template"
            print d_ref_choice, "Reference"
            print name_lig, "Lig Name"
            print name_pose, "Pose Name"
            print "-".join(l_res), "List of residues"
            print "**********************"

        conv = OBConversion()
        conv.SetInFormat("mol2")

        st_ligdock = OBMol()
        st_protdock = OBMol()

        conv.ReadFile(st_ligdock, p_pose)
        conv.ReadFile(st_protdock, p_pocket)

        d_check_ring = getringdict(st_protdock)

        ringinteraction(d_check_choice, d_check_ring, l_res, st_protdock, st_ligdock)
        otherinteractions(d_check_choice, l_res, st_protdock, st_ligdock, l_serial_backbone_test)
        hbonddockprot(d_check_choice, l_res, st_protdock, st_ligdock)

        # complet d_out
        if not str(name_lig) + "--" + str(name_pose) in d_out.keys():
            d_out[str(name_lig) + "--" + str(name_pose)] = {}

        d_out[str(name_lig) + "--" + str(name_pose)]["FPI"] = " ".join(
            [collectbit(d_check_choice, [res]) for res in l_res])
        d_out[str(name_lig) + "--" + str(name_pose)]["score"] = gettcfromdict(d_ref_choice, d_check_choice, l_res)

        if debug == 1:
            print "Out dictionnary", d_out[str(name_lig) + "--" + str(name_pose)]

        del st_protdock
        del st_ligdock

        # write filout
        filout.write(str(name_lig_pose) + "\t" + str(d_out[str(name_lig) + "--" + str(name_pose)]["score"]) + "\t")
        filout.write(" ".join([collectbit(d_check_choice, [res]) for res in l_res]) + "\n")

    filout.close()

    y = time()
    print 'Total time taken %.3f s.' % (y - x)

    return FingerPrintInteraction(p_pocket_ref, p_lig_ref, pr_pose_and_pocket, l_res, p_filout)

# FingerPrintInteraction ("/home/borrel/GPCR_projects/result/No_waters/FPI/SUV/SUV_1_pocket.mol2", "/home/borrel/GPCR_projects/result/No_waters/FPI/SUV/SUV_1.mol2", "/home/borrel/GPCR_projects/result/No_waters/FPI/SUV/Pose/", ['PRO131', 'GLN134', 'ALA110', 'THR111', 'VAL114', 'TRP120', 'ASN324', 'ILE320', 'VAL138', 'ILE130', 'THR135', 'GLN187', 'GLU212', 'VAL353', 'HIS350', 'TYR317', 'PHE227', 'TYR354', 'HIS224'], "/home/borrel/GPCR_projects/result/No_waters/FPI/SUV/FPI_SUV")
"""