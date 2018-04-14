#!/usr/bin/env python
# -*- coding:utf-8 -*-
from time import time, sleep
import getopt, sys, os, copy, glob, re
#from openbabel import OBMol, OBConversion
from optparse import OptionParser
#from bitarray import bitarray
from os import path

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
    #print 'Total time taken %.3f s.' % (y - x)

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

