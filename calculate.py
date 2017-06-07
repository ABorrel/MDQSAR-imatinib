from math import sqrt

def jaccardIndex(lbit1, lbit2):

    if len(lbit1) != len(lbit2):
        print "ERROR: no same length of bits - l.6 calculate"
        return "ERROR"

    i = 0
    identic = 0.0
    while i < len(lbit1):
        if lbit1[i] == lbit2[i]:
            identic += 1.0
        i += 1

    if identic != 0.0:
        score = identic / (len(lbit1) + len(lbit2) - identic)
    else:
        score = 0.0


    return score



def RMSDTwoList(l_atom1, l_atom2):
    nb_ca = 0.0
    d_max = {"value": 0.0}
    diff_position_all = 0.0
    diff_position_ca = 0.0

    if len(l_atom1) != len(l_atom2) or len(l_atom2) == 0:
        print "ERROR - RMSD: list length different or null"
        return []
    else:
        i = 0
        while i < len(l_atom1):
            if l_atom1[i].name != l_atom2[i].name and l_atom1[i].resName != l_atom2[i].resName:
                print l_atom1[i].name, l_atom2[i].name
                print "ERROR"
                return []
            else:

                dist_atoms = l_atom1[i].euclidiendist(l_atom2[i])
                diff_position_all = diff_position_all + dist_atoms

                if l_atom1[i].name == "CA":
                    diff_position_ca = diff_position_ca + dist_atoms
                    nb_ca = nb_ca + 1

                if dist_atoms > d_max["value"]:
                    d_max["value"] = dist_atoms
                    d_max["atom"] = l_atom1[i].name + "-" + l_atom2[i].name + "_" + l_atom1[i].resName + "-" + \
                                    l_atom2[i].resName

            i = i + 1
            #     print d_max

    RMSDall = sqrt(diff_position_all / len(l_atom1))
    if nb_ca >0:
        RMSDCA = sqrt(diff_position_ca / nb_ca)
    else:
        RMSDCA = "NA"

    return [RMSDall, RMSDCA, d_max["value"], len(l_atom1)]

#print jaccardIndex("001011101", "001001001")



