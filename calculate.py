

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



#print jaccardIndex("001011101", "001001001")



