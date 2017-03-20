


def mergeDict(l_dict):

    dout = {}

    lk = l_dict[0].keys()

    for k in lk:
        lval = [l_dict[i][k] for i in range(0,len(l_dict))]
        #identic
        if lval.count(lval[0]) == len(lval):
            dout[k] = lval[0]
        else:
            dout[k] = "----".join(lval)

    return dout




def convertUnit(l_values, l_units):
    """Convert list of affinity un uM"""

    #print l_values
    #print l_units
    lout = []
    i = 0
    imax = len(l_values)
    while i < imax:
        if l_units[i] == "uM" or l_units[i] == "10'-6M":
            lout.append(l_values[i])
            i += 1
        elif l_units[i] == "nM":
            val = float(l_values[i])
            val = val/1000
            lout.append(val)
            i += 1
        else:
            ffff

    #print lout

    return lout