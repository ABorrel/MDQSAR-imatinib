# -*- coding: utf-8 -*-
"""
##############################################################################
The calculation of whole holistic invariant molecular descriptors (WHIM).
You can get 70 molecular decriptors. You can freely use and distribute it.
If you hava any problem, you could contact with us timely!
Authors: Dongsheng Cao and Yizeng Liang, Qingsong Xu
Date: 2011.04.19
Email: oriental-cds@163.com
##############################################################################
"""
import scipy
import scipy.linalg
from AtomProperty import GetRelativeAtomicProperty



Version=1.0
#############################################################################




def XPreCenter(X):
    """
    #################################################################
    Center the data matrix X
    #################################################################
    """
    Xdim=scipy.size(X,axis=0)
    Xmean=scipy.mean(X,axis=0)
    Xmean=scipy.matrix(Xmean)
    Xp=X-scipy.ones([Xdim,1])*Xmean
   
    return Xp


def GetPropertyMatrix(AtomLabel,proname='m'):
    """
    #################################################################
    #################################################################
    """
    res=[]
    for i in AtomLabel:
        res.append(GetRelativeAtomicProperty(i,proname))
    
    return scipy.matrix(scipy.diag(res))
        


def GetSVDEig(CoordinateMatrix,AtomLabel,proname='u'):
    """
    """
    nAtom,kc= CoordinateMatrix.shape   
    
    if proname=='u':
        weight=scipy.matrix(scipy.eye(nAtom))
    else:      
        weight=GetPropertyMatrix(AtomLabel,proname)
    
    S=XPreCenter(CoordinateMatrix)

    u,s,v=scipy.linalg.svd(S.T*weight*S/sum(scipy.diag(weight)))
    
    return s




def GetWHIM1(CoordinateMatrix,AtomLabel,proname='u'):
    """
    #################################################################
    WHIM descriptors
    --->L1u
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    return round(s[0],3)


def GetWHIM2(CoordinateMatrix,AtomLabel,proname='u'):
    """
    #################################################################
    WHIM descriptors
    
    --->L2u
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    return round(s[1],3)


def GetWHIM3(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->L3u
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    return round(s[2],3)


def GetWHIM4(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->Tu
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    T=round(sum(s),3)
    
    return T

def GetWHIM5(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->Au
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    A=s[0]*s[1]+s[0]*s[2]+s[1]*s[2]
    
    return round(A,3)


def GetWHIM6(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->Vu
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    A=s[0]*s[1]+s[0]*s[2]+s[1]*s[2]
    T=sum(s)
    V=A+T+s[0]*s[1]*s[2]
    
    return round(V,3)

def GetWHIM7(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->P1u
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    return round(s[0]/(s[0]+s[1]+s[2]),3)


def GetWHIM8(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->P2u
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    return round(s[1]/(s[0]+s[1]+s[2]),3)




def GetWHIM9(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->Ku
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    res=0.0
    for i in s:
        res=res+abs(i/sum(s)-1/3.0)
        
    Ku=3.0/4*res
    
    return round(Ku,3) 


def GetWHIM10(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->E1u
    #################################################################
    """
    
    nAtom,kc= CoordinateMatrix.shape   
    
    if proname=='u':
        weight=scipy.matrix(scipy.eye(nAtom))
    else:      
        weight=GetPropertyMatrix(AtomLabel,proname)
    S=XPreCenter(CoordinateMatrix)
    u,s,v=scipy.linalg.svd(S.T*weight*S/sum(scipy.diag(weight)))
     
    res=scipy.power(s[0],2)*nAtom/sum(scipy.power(S*scipy.matrix(u[:,0]).T,4))
    
    return round(float(res.real),3)
    
    
def GetWHIM11(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->E2u
    #################################################################
    """
    nAtom,kc= CoordinateMatrix.shape   
    
    if proname=='u':
        weight=scipy.matrix(scipy.eye(nAtom))
    else:      
        weight=GetPropertyMatrix(AtomLabel,proname)
    S=XPreCenter(CoordinateMatrix)
    u,s,v=scipy.linalg.svd(S.T*weight*S/sum(scipy.diag(weight)))
     
    res=scipy.power(s[1],2)*nAtom/sum(scipy.power(S*scipy.matrix(u[:,1]).T,4))
    
    return round(float(res.real),3)
    

def GetWHIM12(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->E3u
    #################################################################
    """
    nAtom,kc= CoordinateMatrix.shape   
    
    if proname=='u':
        weight=scipy.matrix(scipy.eye(nAtom))
    else:      
        weight=GetPropertyMatrix(AtomLabel,proname)
    S=XPreCenter(CoordinateMatrix)
    u,s,v=scipy.linalg.svd(S.T*weight*S/sum(scipy.diag(weight)))
     
    try:
        res=scipy.power(s[2],2)*nAtom/sum(scipy.power(S*scipy.matrix(u[:,2]).T,4))
        return round(float(res.real),3)
    except: return "NA"
    
def GetWHIM13(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->Du
    #################################################################
    """
    c1=GetWHIM10(CoordinateMatrix,AtomLabel,proname)
    c2=GetWHIM11(CoordinateMatrix,AtomLabel,proname)
    c3=GetWHIM12(CoordinateMatrix,AtomLabel,proname)
    Du=c1+c2+c3
    
    return round(float(Du),3)




def GetWHIM14(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors

    --->P3u
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)

    return round(s[2]/(s[0]+s[1]+s[2]),3)

