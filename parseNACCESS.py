'''
Created on 13 sept. 2013

@author: Borrel Alexandre
'''


from re import search


def fileRSA(path_filin_rsa):
    
    list_res_out = []
    filin = open (path_filin_rsa, "r")
    list_lines_rsa = filin.readlines()
    filin.close()
    
    for line_rsa in list_lines_rsa : 
        if search("^RES", line_rsa) :
            list_res_out.append(parseLineRsa(line_rsa))
    
    return list_res_out


def fileASA(pasa):
    """
    Parse global file asa
    args: -> path file asa
    return: file parsed
    """
    
    
    list_atom_out = []
    filin = open(pasa, "r")
    list_lines_rsa = filin.readlines()
    filin.close()
    
    for line_rsa in list_lines_rsa : 
        if search("^ATOM", line_rsa) :
            list_atom_out.append (parseLineAsa(line_rsa))
    return list_atom_out


def parseLineRsa(lineRes):
    """
    Parse line NACCESS, retrieve residue ID and solvent exposition
    args: line in file
    out: dictionary with residue ID and solvent exposition
    """
    
    out = {}
    
    out["resSeq"] = int(lineRes[9:13].replace(" ", ""))
    out["ABS"] = float(lineRes[16:22].replace(" ", ""))
    out["REL"] = float(lineRes[23:28].replace(" ", ""))
    out["chainID"] = str(lineRes[8])
    out["line"] = str(lineRes)
    return out


def parseLineAsa(lineAtom):
    """
    Parse line NACCESS, retrieve residue ID and solvent exposition
    args: -> line in file
    out: -> dictionary parsed line
    """
    
    out = {}
    out["atomSeq"] = int(lineAtom[4:11].replace(" ",""))
    out["name"] = lineAtom[12:16].replace (" ", "")
    out["resName"] = lineAtom[12:16].replace (" ", "")
    out["chainID"] = str(lineAtom[21])
    out["ABS"] = float(lineAtom[56:63].replace(" ",""))
    out["resSeq"] = int (lineAtom[22:26].replace (" ", ""))
    out["radius"] = float (lineAtom[64:70].replace (" ",""))
    out["line"] = str(lineAtom)
    return out
