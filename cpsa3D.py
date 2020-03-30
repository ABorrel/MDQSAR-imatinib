import asa
import descriptors3D


def GetChargeSA(lcoordinates, RadiusProbe=1.5, n_sphere_point=960):
    """Get the list form for all atoms in a molecule.
    It includes the atom symbol, charge and partial solvent-accessible
    surface areas.
    Note that this is list form whose element is still list form of each atom.

    """
    atoms=descriptors3D.GetAtomClassList(lcoordinates)
    FASA=asa.calculate_asa(atoms, RadiusProbe, n_sphere_point)
    res=[]
    for i in range(len(FASA)):
        res.append([lcoordinates[i][3],lcoordinates[i][4],FASA[i]])
    return res


def CalculateASA(ChargeSA):
    """ The calculation of solvent-accessible surface areas
    -->ASA """
    res=0.0
    for i in ChargeSA:
        res=res+i[2]
    return res


def CalculateMSA(lcoordinates):
    """ The calculation of molecular surface areas
    -->MSA"""
    ChargeSA=GetChargeSA(lcoordinates, RadiusProbe=0, n_sphere_point=960)
    res=0.0
    for i in ChargeSA:
        res=res+i[2]
    return res


def CalculatePNSA1(ChargeSA):
    """ The calculation of partial negative area
    It is the sum of the solvent-accessible surface areas of all 
    negatively charged atoms.
    -->PNSA1"""
    res=0.0
    for i in ChargeSA:
        if float(i[1])<0:
            res=res+i[2]

    return res


def CalculatePPSA1(ChargeSA):
    """ The calculation of partial negative area
    It is the sum of the solvent-accessible surface areas of
    all positively charged atoms.
    -->PPSA1
    """
    res=0.0
    for i in ChargeSA:
        if float(i[1])>0:
            res=res+i[2]

    return res

def CalculatePNSA2(ChargeSA):
    """The calculation of total charge wighted negative surface area
    It is the partial negative solvent-accessible surface area
    multiplied by the total negative charge.
    -->PNSA2
    """
    temp1=temp2=0.0
    for i in ChargeSA:
        if float(i[1])<0:
            temp1=temp1+float(i[1])
            temp2=temp2+i[2]
    res=temp1*temp2

    return res


def CalculatePPSA2(ChargeSA):
    """The calculation of total charge wighted negative surface area
    It is the partial negative solvent-accessible surface area 
    multiplied by the total positive charge.
    -->PPSA2 """

    temp1=temp2=0.0
    for i in ChargeSA:
        if float(i[1])>0:
            temp1=temp1+float(i[1])
            temp2=temp2+i[2]
    res=temp1*temp2

    return res

def CalculatePNSA3(ChargeSA):
    """The calculation of atom charge weighted negative surface ares
    It is the sum of the products of atomic solvent-accessible 
    surface area and partial charges over all negatively charges atoms.
    -->PNSA3 """

    res=0.0
    for i in ChargeSA:
        if float(i[1])<0:
            res=res+float(i[1])*i[2]
    return res


def CalculatePPSA3(ChargeSA):
    """The calculation of atom charge weighted positive surface ares
    It is the sum of the products of atomic solvent-accessible
    surface area and partial charges over all positively charges atoms.
    -->PPSA3"""

    res=0.0
    for i in ChargeSA:
        if float(i[1])>0:
            res=res+float(i[1])*i[2]
    return res


def CalculateDPSA1(ChargeSA):
    """ The calculation of difference in charged partial surface area
    -->DPSA1"""
    return CalculatePPSA1(ChargeSA)-CalculatePNSA1(ChargeSA)


def CalculateDPSA2(ChargeSA):
    """The calculation of difference in total charge weighted partial
    surface area
    -->DPSA2"""
    return CalculatePPSA2(ChargeSA)-CalculatePNSA2(ChargeSA)


def CalculateDPSA3(ChargeSA):
    """The calculation of difference in atomic charge weighted surface area
    -->DPSA3"""
    return CalculatePPSA3(ChargeSA)-CalculatePNSA3(ChargeSA)

def CalculateFNSA1(ChargeSA):
    """The calculation of fractional charged partial negative surface areas
    -->FNSA1"""
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]

    return CalculatePNSA1(ChargeSA)/temp


def CalculateFNSA2(ChargeSA):
    """The calculation of fractional charged partial negative surface areas
    -->FNSA2
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]

    return CalculatePNSA2(ChargeSA)/temp


def CalculateFNSA3(ChargeSA):
    """The calculation of fractional charged partial negative surface areas
    -->FNSA3
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]

    return CalculatePNSA3(ChargeSA)/temp


def CalculateFPSA1(ChargeSA):
    """The calculation of fractional charged partial negative surface areas
    -->FPSA1
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]

    return CalculatePPSA1(ChargeSA)/temp


def CalculateFPSA2(ChargeSA):
    """The calculation of fractional charged partial negative surface areas
    -->FPSA2
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]

    return CalculatePPSA2(ChargeSA)/temp


def CalculateFPSA3(ChargeSA):
    """The calculation of fractional charged partial negative surface
    areas
    -->FPSA3
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]

    return CalculatePPSA3(ChargeSA)/temp


def CalculateWNSA1(ChargeSA):
    """The calculation of surface weighted charged partial negative 
    surface areas
    -->WNSA1
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]

    return CalculatePNSA1(ChargeSA)*temp/1000


def CalculateWNSA2(ChargeSA):
    """ The calculation of surface weighted charged partial negative 
    surface areas
    -->WNSA2
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]

    return CalculatePNSA2(ChargeSA)*temp/1000


def CalculateWNSA3(ChargeSA):
    """The calculation of surface weighted charged partial negative 
    surface areas
    -->WNSA3
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]

    return CalculatePNSA3(ChargeSA)*temp/1000


def CalculateWPSA1(ChargeSA):
    """The calculation of surface weighted charged partial negative 
    surface areas
    -->WPSA1
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]

    return CalculatePPSA1(ChargeSA)*temp/1000


def CalculateWPSA2(ChargeSA):
    """The calculation of surface weighted charged partial negative
    surface areas
    -->WPSA2
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]

    return CalculatePPSA2(ChargeSA)*temp/1000


def CalculateWPSA3(ChargeSA):
    """The calculation of surface weighted charged partial negative
    surface areas
    -->WPSA3
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]

    return CalculatePPSA3(ChargeSA)*temp/1000


def CalculateTASA(ChargeSA):
    """The calculation of total hydrophobic surface area
    -->TASA
    """
    res=0.0
    for i in ChargeSA:
        if abs(float(i[1]))<0.2:
            res=res+i[2]
    return res


def CalculateTPSA(ChargeSA):
    """The calculation of total polar surface area
    -->PSA
    """
    res=0.0
    for i in ChargeSA:
        if abs(float(i[1]))>=0.2:
            res=res+i[2]
    return res


def CalculateFractionTATP(ChargeSA):
    """The fraction between TASA and TPSA
    --->FrTATP
    """
    res=0.0
    if CalculateTPSA(ChargeSA)==0:
        return res
    else:
        return CalculateTASA(ChargeSA)/CalculateTPSA(ChargeSA)


def CalculateRASA(ChargeSA):
    """The calculation of relative hydrophobic surface area
    -->RASA
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
    return CalculateTASA(ChargeSA)/temp


def CalculateRPSA(ChargeSA):
    """The calculation of relative polar surface area
    -->RPSA
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
    return CalculateTPSA(ChargeSA)/temp


def CalculateRNCS(ChargeSA):
    """The calculation of relative negative charge surface area
    -->RNCS
    """
    charge=[]
    for i in ChargeSA:
        charge.append(float(i[1]))

    temp=[]
    for i in ChargeSA:
        temp.append(i[2])

    try:
        RNCG=min(charge)/sum([i for i in charge if i<0])
        return  temp[charge.index(min(charge))]/RNCG
    except:
        return "NA"


def CalculateRPCS(ChargeSA):
    """The calculation of relative positive charge surface area
    -->RPCS
    """
    charge=[]
    for i in ChargeSA:
        charge.append(float(i[1]))

    temp=[]
    for i in ChargeSA:
        temp.append(i[2])

    try:
        RPCG=max(charge)/sum([i for i in charge if i>0])
        return  temp[charge.index(min(charge))]/RPCG
    except:
        return "NA"

