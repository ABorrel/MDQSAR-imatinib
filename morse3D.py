# -*- coding: utf-8 -*-
"""
##############################################################################

The calculation of 3D morse descriptors. You can get 210 molecular

decriptors. You can freely use and distribute it. If you hava  

any problem, you could contact with us timely!

Authors: Dongsheng Cao and Yizeng Liang.

Date: 2012.11.13

Email: oriental-cds@163.com

##############################################################################
"""

from AtomProperty import GetRelativeAtomicProperty
import scipy
import math
import descriptors3D

_beta=100


def _GetR(n=32):
    """
    *Internal Use Only*
    Obtain the parameter R in RDF equation.
    """
    R=[]
    for i in range(1,n+1):
        R.append(float(i*1))
    return R


def _GetAtomDistance(x,y):
    """
    *Internal Use Only*
    Obtain the Elucidian distance based on the coordinates of two atoms
    """

    temp=[math.pow(x[0]-y[0],2),math.pow(x[1]-y[1],2),math.pow(x[2]-y[2],2)]
    res=math.sqrt(sum(temp))
    return res
    

def _GetGementricalDistanceMatrix(CoordinateList):
    """
    *Internal Use Only*
    Obtain the distance matrix of a molecule based on coordinate list
    """
    NAtom=len(CoordinateList)
    DistanceMatrix=scipy.zeros((NAtom,NAtom))
    for i in range(NAtom-1):
        for j in range(i+1,NAtom):
            DistanceMatrix[i,j]=_GetAtomDistance(CoordinateList[i],CoordinateList[j])
            DistanceMatrix[j,i]=DistanceMatrix[i,j]
    return DistanceMatrix

    
def CalculateUnweightMoRSE(lcoordinates):
    """
    The calculation of  unweighted 3-D MoRse descriptors
    """
    R=_GetR(n=30)
    temp=[]
#    ChargeCoordinates=_ReadCoordinates('temp.arc')
    for i in lcoordinates:
        #if i[0]!='H':
        temp.append([float(i[0]),float(i[1]),float(i[2])])
        
    DM=_GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}
    
    for kkk,Ri in enumerate(R):        
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+math.sin(Ri*DM[j,k])/(Ri*DM[j,k])
        RDFresult['MoRSE'+'U'+str(kkk+1)]=round(res,3)
        
    return RDFresult


def CalculateChargeMoRSE(lcoordinates):
    
    """
    The calculation of  3-D MoRse descriptors
    based on atomic charge.
    """
    R=_GetR(n=30)
    temp=[]
    charge=[]
    for i in lcoordinates:
        #if i[0]!='H':
        temp.append([float(i[0]),float(i[1]),float(i[2])])
        charge.append(float(i[4]))
    DM=_GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}
    
    for kkk,Ri in enumerate(R):        
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+charge[j]*charge[k]*math.sin(Ri*DM[j,k])/(Ri*DM[j,k])
        RDFresult['MoRSE'+'C'+str(kkk+1)]=round(res,3)
        
    return RDFresult



def CalculateMassMoRSE(lcoordinantes):
    """
    #################################################################
    The calculation of  3-D MoRse descriptors 
    
    based on atomic mass.
    #################################################################
    """
    R=_GetR(n=30)
    temp=[]
    mass=[descriptors3D.get_atomicMass(coords[3]) for coords in lcoordinantes]
    for i in lcoordinantes:
        #if i[0]!='H':
        temp.append([float(i[0]),float(i[1]),float(i[2])])
    DM=_GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}
    
    for kkk,Ri in enumerate(R):        
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+mass[j]*mass[k]*math.sin(Ri*DM[j,k])/(Ri*DM[j,k])
        RDFresult['MoRSE'+'M'+str(kkk+1)]=round(res/144,3)
        
    return RDFresult    



def CalculateAtomicNumberMoRSE(lcoordinates):
    """
    #################################################################
    The calculation of  3-D MoRse descriptors 
    
    based on atomic number.
    #################################################################
    """
    R=_GetR(n=30)
    temp=[]
    mass = [descriptors3D.get_atomicMass(coords[3]) for coords in lcoordinates]
    for i in lcoordinates:
        #if i[0]!='H':
        temp.append([float(i[0]),float(i[1]),float(i[2])])
    DM=_GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}
    
    for kkk,Ri in enumerate(R):        
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+mass[j]*mass[k]*math.sin(Ri*DM[j,k])/(Ri*DM[j,k])
        RDFresult['MoRSE'+'N'+str(kkk+1)]=round(res/144,3)
        
    return RDFresult       




def CalculatePolarizabilityMoRSE(lcoordinates):
    """
    #################################################################
    The calculation of  3-D MoRse descriptors 
    
    based on atomic polarizablity.
    #################################################################
    """
    R=_GetR(n=30)
    temp=[]
    polarizability=[]
    for i in lcoordinates:
        #if i[0]!='H':
        temp.append([float(i[0]),float(i[1]),float(i[2])])
        polarizability.append(GetRelativeAtomicProperty(i[3],'alapha'))
    DM=_GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}
    
    for kkk,Ri in enumerate(R):        
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+polarizability[j]*polarizability[k]*math.sin(Ri*DM[j,k])/(Ri*DM[j,k])
        RDFresult['MoRSE'+'P'+str(kkk+1)]=round(res,3)
        
    return RDFresult



def CalculateSandersonElectronegativityMoRSE(lcoordinates):
    """
    #################################################################
    The calculation of  3-D MoRse descriptors 
    
    based on atomic sanderson electronegativity.
    #################################################################
    """
    R=_GetR(n=30)
    temp=[]
    En=[]
    for i in lcoordinates:
        #if i[0]!='H':
        temp.append([float(i[0]),float(i[1]),float(i[2])])
        En.append(GetRelativeAtomicProperty(i[3],'En'))
    DM=_GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}
    
    for kkk,Ri in enumerate(R):        
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+En[j]*En[k]*math.sin(Ri*DM[j,k])/(Ri*DM[j,k])
        RDFresult['MoRSE'+'E'+str(kkk+1)]=round(res,3)
        
    return RDFresult



def CalculateVDWVolumeMoRSE(ChargeCoordinates):
    """
    #################################################################
    The calculation of  3-D MoRse descriptors 
    
    based on atomic van der Waals volume.
    #################################################################
    """
    R=_GetR(n=30)
    temp=[]
    VDW=[]
    for i in ChargeCoordinates:
        #if i[0]!='H':
        temp.append([float(i[0]),float(i[1]),float(i[2])])
        VDW.append(GetRelativeAtomicProperty(i[3],'V'))
    DM=_GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}
    
    for kkk,Ri in enumerate(R):        
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+VDW[j]*VDW[k]*math.sin(Ri*DM[j,k])/(Ri*DM[j,k])
        RDFresult['MoRSE'+'V'+str(kkk+1)]=round(res,3)
        
    return RDFresult

