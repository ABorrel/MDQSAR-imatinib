from AtomProperty import GetRelativeAtomicProperty
import scipy
import math
import descriptors3D


#########################################################################
###set the parameters in RDF equation
_beta=100
#########################################################################

def _GetR(n=30):
    """Obtain the parameter R in RDF equation.
    """
    R=[]
    for i in range(2,n+2):
        R.append(float(i*0.5))
    return R


def GetAtomDistance(x,y):
    """
    Obtain the Elucidian distance based on the
    coordinates of two atoms
    """

    temp=[math.pow(x[0]-y[0],2),math.pow(x[1]-y[1],2),math.pow(x[2]-y[2],2)]
    res=math.sqrt(sum(temp))
    return res
    

def GetGementricalDistanceMatrix(CoordinateList):
    """
    Obtain the distance matrix of a molecule based
    on coordinate list
    """
    NAtom=len(CoordinateList)
    DistanceMatrix=scipy.zeros((NAtom,NAtom))
    for i in range(NAtom-1):
        for j in range(i+1,NAtom):
            DistanceMatrix[i,j]=GetAtomDistance(CoordinateList[i],CoordinateList[j])
            DistanceMatrix[j,i]=DistanceMatrix[i,j]
    return DistanceMatrix

    
def CalculateUnweightRDF(lcoordinates):
    """
    The calculation of unweighted radial distribution
    function (RDF) descriptors.
    """
    R=_GetR(n=30)
    temp=[]

    for i in lcoordinates:
        #if i[0]!='H': Have to considerate the H?
        temp.append([float(i[0]),float(i[1]),float(i[2])])

    DM=GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}

    for kkk,Ri in enumerate(R):
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+math.exp(-_beta*math.pow(Ri-DM[j,k],2))
        RDFresult['RDF'+'U'+str(kkk+1)]=round(res,3)

    return RDFresult



def CalculateChargeRDF(lcoordinates):
    """
    The calculation of  radial distribution function
    (RDF) descriptors based on atomic charge.
    """
    R=_GetR(n=30)
    temp=[]
    Charge=[]
    for i in lcoordinates:
        #if i[0]!='H':
        temp.append([float(i[0]),float(i[1]),float(i[2])])
        Charge.append(float(i[4]))

    DM=GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}

    for kkk,Ri in enumerate(R):
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+Charge[j]*Charge[k]*math.exp(-_beta*math.pow(Ri-DM[j,k],2))
        RDFresult['RDF'+'C'+str(kkk+1)]=round(res,3)

    return RDFresult


def CalculateMassRDF(lcoordinates):
    """
    The calculation of radial distribution function (RDF)
    descriptors based on atomic mass.
    """
    mass=[descriptors3D.get_atomicMass(coords[3]) for coords in lcoordinates]
    R=_GetR(n=30)
    temp=[]
    for i in lcoordinates:
        #if i[0]!='H':
        temp.append([float(i[0]),float(i[1]),float(i[2])])

    DM=GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}

    for kkk,Ri in enumerate(R):
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+mass[j]*mass[k]*math.exp(-_beta*math.pow(Ri-DM[j,k],2))
        RDFresult['RDF'+'M'+str(kkk+1)]=round(res/144,3)

    return RDFresult



def CalculatePolarizabilityRDF(lcoordinates):
    """
    The calculation of  radial distribution function
    (RDF) descriptors based on atomic polarizability.
    """
    R=_GetR(n=30)
    temp=[]
    polarizability=[]
#    lcoordinates=_ReadCoordinates('temp.arc')
    for i in lcoordinates:
        #if i[0]!='H':
        temp.append([float(i[0]),float(i[1]),float(i[2])])
        polarizability.append(GetRelativeAtomicProperty(i[3],'alapha'))

    DM=GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}

    for kkk,Ri in enumerate(R):
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+polarizability[j]*polarizability[k]*math.exp(-_beta*math.pow(Ri-DM[j,k],2))
        RDFresult['RDF'+'P'+str(kkk+1)]=round(res,3)

    return RDFresult



def CalculateSandersonElectronegativityRDF(lcoordinates):
    """
    The calculation of  radial distribution function
    (RDF) descriptors based on atomic electronegativity.
    """
    R=_GetR(n=30)
    temp=[]
    EN=[]
#    lcoordinates=_ReadCoordinates('temp.arc')
    for i in lcoordinates:
        #if i[0]!='H':
        temp.append([float(i[0]),float(i[1]),float(i[2])])
        EN.append(GetRelativeAtomicProperty(i[3],'En'))

    DM=GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}

    for kkk,Ri in enumerate(R):
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+EN[j]*EN[k]*math.exp(-_beta*math.pow(Ri-DM[j,k],2))
        RDFresult['RDF'+'E'+str(kkk+1)]=round(res,3)

    return RDFresult


def CalculateVDWVolumeRDF(lcoordinates):
    """
    The calculation of  radial distribution function
    (RDF) descriptors based on atomic van der Waals volume.
    """
    R=_GetR(n=30)
    temp=[]
    VDW=[]
#    lcoordinates=_ReadCoordinates('temp.arc')
    for i in lcoordinates:
        #if i[0]!='H':
        temp.append([float(i[0]),float(i[1]),float(i[2])])
        VDW.append(GetRelativeAtomicProperty(i[3],'V'))

    DM=GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}

    for kkk,Ri in enumerate(R):
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+VDW[j]*VDW[k]*math.exp(-_beta*math.pow(Ri-DM[j,k],2))
        RDFresult['RDF'+'V'+str(kkk+1)]=round(res,3)

    return RDFresult



def GetRDFUnweighed(mol):
    """
    Obtain all Unweighed radial distribution function descriptors.
    """

    filename='temp'
    lcoordinates=_ReadCoordinates(filename) 
    result=CalculateUnweightRDF(lcoordinates)
     
    return result



def GetRDFCharge(mol):
    """
    Obtain all radial distribution function descriptors based
    on Charge schems.
    """

    filename='temp'
    lcoordinates=_ReadCoordinates(filename) 
    result=CalculateChargeRDF(lcoordinates)
     
    return result
    
    
def GetRDFMass(mol):
    """
    Obtain all radial distribution function descriptors based
    on Mass schems.
    """

    filename='temp'
    lcoordinates=_ReadCoordinates(filename) 
    result=CalculateMassRDF(mol,lcoordinates)
     
    return result
    
    
def GetRDFPolarizability(mol):
    """
    Obtain all radial distribution function descriptors based
    on Polarizability schems.
    """

    filename='temp'
    lcoordinates=_ReadCoordinates(filename) 
    result=CalculatePolarizabilityRDF(lcoordinates)
     
    return result



def GetRDFSandersonElectronegativity(mol):
    """
    Obtain all radial distribution function descriptors based
    onSanderson Electronegativity schems.
    """

    filename='temp'
    lcoordinates=_ReadCoordinates(filename) 
    result=CalculateSandersonElectronegativityRDF(lcoordinates)
     
    return result


def GetRDFVDWVolume(mol):
    """
    Obtain all radial distribution function descriptors based
    on VDW Volume schems.
    """

    filename='temp'
    lcoordinates=_ReadCoordinates(filename) 
    result=CalculateVDWVolumeRDF(lcoordinates)
     
    return result



