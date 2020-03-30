import geo3D
import cpsa3D
import rdf3D
import morse3D
import vector3d
import whim3D

import scipy
from re import search
from os import path
# compute from CHEMPY
# replace MOPAC file by sdf parsing


l3D = ['RDFC6', 'MoRSEN11', 'RDFU8', 'RDFU9', 'RDFU2', 'RDFU3', 'MoRSEN5', 'RDFU1', 'RDFU6', 'RDFU7', 'RDFU4', 'RDFU5',
       'Harary3D', 'P2u', 'MoRSEM6', 'MoRSEM7', 'MoRSEM4', 'MoRSEM5', 'MoRSEM2', 'MoRSEM3', 'MoRSEE30', 'MoRSEM1',
       'MoRSEN4', 'MoRSEM8', 'MoRSEM9', 'MoRSEU10', 'MoRSEU11', 'MoRSEU12', 'MoRSEU13', 'MoRSEU14', 'MoRSEU15',
       'MoRSEU16', 'MoRSEU17', 'MoRSEU18', 'MoRSEU19', 'FPSA3', 'FPSA2', 'FPSA1', 'GeDi', 'MoRSEV19', 'MoRSEN10',
       'MoRSEV13', 'SPAN', 'MoRSEV11', 'MoRSEV10', 'MoRSEV17', 'MoRSEV16', 'MoRSEV15', 'MoRSEN16', 'RDFM14', 'RDFM15',
       'RDFM16', 'RDFE9', 'RDFM10', 'RDFM11', 'RDFM12', 'RASA', 'RDFE2', 'RDFE3', 'MoRSEC25', 'MoRSEC24', 'RDFM18',
       'RDFM19', 'grav', 'RDFE5', 'WNSA1', 'WNSA2', 'WNSA3', 'L2p', 'RDFP15', 'RDFP14', 'RDFP17', 'RDFP16', 'RDFP11',
       'RDFP10', 'RDFP13', 'RDFP12', 'MoRSEP30', 'RDFP19', 'RDFP18', 'E2p', 'Dm', 'P3e', 'MoRSEM18', 'MoRSEM19',
       'Petitj3D', 'MoRSEM10', 'MoRSEM11', 'MoRSEM12', 'MoRSEM13', 'MoRSEM14', 'MoRSEM15', 'MoRSEM16', 'MoRSEM17',
       'RDFC27', 'RDFC26', 'RDFC25', 'RDFC24', 'RDFC23', 'RDFC22', 'RDFC21', 'RDFC20', 'MoRSEC28', 'RDFU30', 'RDFC29',
       'RDFC28', 'MoRSEU30', 'L1u', 'L1v', 'L2v', 'L1p', 'RDFP5', 'RDFP4', 'RDFP7', 'RDFP6', 'RDFP1', 'RDFP3', 'RDFP2',
       'L1e', 'RDFP9', 'RDFP8', 'MoRSEP5', 'P1e', 'MoRSEP4', 'PSA', 'MoRSEP7', 'P1p', 'MoRSEP6', 'RDFE18', 'RDFE19',
       'P1v', 'RDFE14', 'RDFE15', 'RDFE16', 'RDFE17', 'PPSA1', 'RDFE11', 'RDFE12', 'PPSA2', 'MoRSEP11', 'MoRSEP10',
       'MoRSEP13', 'MoRSEP12', 'RPCS', 'MoRSEP14', 'MoRSEN9', 'MoRSEN8', 'DPSA1', 'MoRSEC30', 'DPSA3', 'DPSA2',
       'MoRSEN3', 'MoRSEN2', 'MoRSEN1', 'RDFP30', 'E2e', 'MoRSEN17', 'L3e', 'TASA', 'RDFC19', 'MoRSEV14', 'MoRSEM30',
       'MoRSEP8', 'L3v', 'RDFC16', 'L3u', 'RDFV30', 'L3p', 'RDFC14', 'W3DH', 'RDFC15', 'MoRSEC23', 'MoRSEN15',
       'MoRSEP16', 'RPSA', 'P3m', 'MEcc', 'MoRSEC22', 'MoRSEN14', 'MoRSEP1', 'MoRSEN23', 'P3p', 'P3v', 'MoRSEP19',
       'P3u', 'RDFV7', 'RDFC18', 'RDFV6', 'FNSA1', 'RDFC17', 'FNSA3', 'FNSA2', 'RDFC12', 'RDFC13', 'RDFC10', 'RDFC11',
       'P2p', 'RDFV4', 'MoRSEP22', 'RDFV3', 'MoRSEP18', 'RDFV2', 'RDFU21', 'RDFU20', 'RDFU23', 'RDFU22', 'RDFU25',
       'RDFM4', 'RDFU27', 'RDFU26', 'RDFU29', 'RDFU28', 'MoRSEN28', 'MoRSEN29', 'RDFV19', 'RDFV18', 'WPSA2', 'RDFV16',
       'RDFV15', 'RDFV14', 'RDFV13', 'RDFV12', 'RDFV11', 'RDFV10', 'MoRSEN26', 'MoRSEP23', 'MoRSEN27', 'MoRSEN24',
       'MoRSEP25', 'MoRSEN25', 'MoRSEP26', 'MoRSEE8', 'MoRSEE9', 'MoRSEE6', 'MoRSEN22', 'MoRSEE4', 'MoRSEP27',
       'MoRSEE2', 'MoRSEE3', 'RDFU24', 'MoRSEN20', 'MoRSEC9', 'ASPAN', 'RDFE10', 'MoRSEN21', 'Te', 'Vm', 'Vp',
       'MoRSEV18', 'PPSA3', 'Vv', 'RDFE13', 'E2u', 'RDFC30', 'E2v', 'P1m', 'MoRSEV12', 'MoRSEP15', 'MoRSEP17',
       'MoRSEU8', 'MoRSEU9', 'MoRSEU6', 'MoRSEU7', 'MoRSEU4', 'MoRSEU5', 'MoRSEU2', 'MoRSEU3', 'MoRSEU1', 'RDFM29',
       'RDFM28', 'MoRSEN6', 'RDFM21', 'RDFM20', 'RDFM23', 'RDFM22', 'RDFM25', 'RDFM24', 'RDFM27', 'RDFM26', 'RDFM2',
       'RDFM3', 'RDFV5', 'RDFM1', 'RDFM6', 'RDFM7', 'RDFV1', 'RDFM5', 'MoRSEP20', 'MoRSEP21', 'RDFM8', 'RDFM9',
       'MoRSEP24', 'RDFE8', 'RDFV9', 'RDFV8', 'MoRSEC8', 'RDFV17', 'RDFM17', 'WPSA3', 'AGDD', 'MoRSEC1', 'MoRSEC2',
       'MoRSEC3', 'MoRSEC4', 'MoRSEC5', 'MoRSEC6', 'MoRSEC7', 'MoRSEE29', 'MoRSEE28', 'WPSA1', 'MoRSEC29', 'MoRSEE21',
       'MoRSEE20', 'MoRSEE23', 'RDFM13', 'MoRSEE25', 'MoRSEE24', 'MoRSEE27', 'MoRSEE26', 'MoRSEC27', 'MoRSEC26', 'Ae',
       'RDFE1', 'RDFE6', 'RDFE7', 'RDFE4', 'MoRSEV28', 'MoRSEV29', 'MoRSEV26', 'MoRSEV27', 'MoRSEV24', 'MoRSEC20',
       'MoRSEV22', 'MoRSEV23', 'MoRSEV20', 'MoRSEV21', 'MoRSEC12', 'MoRSEC13', 'MoRSEC10', 'MoRSEC11', 'MoRSEC16',
       'MoRSEC17', 'MoRSEC14', 'MoRSEC15', 'P2m', 'MoRSEC18', 'MoRSEC19', 'RDFE30', 'RDFE21', 'RDFE20', 'RDFE23',
       'RDFE22', 'RDFE25', 'RDFE24', 'RDFE27', 'RDFE26', 'RDFE29', 'RDFE28', 'RDFP20', 'RDFP21', 'RDFP22', 'RDFP23',
       'RDFP24', 'RDFP25', 'RDFP26', 'RDFP27', 'RDFP28', 'RDFP29', 'MoRSEV6', 'MoRSEN13', 'MoRSEV30', 'Dv', 'RDFV26',
       'RDFV27', 'RDFV24', 'RDFV25', 'RDFV22', 'RDFV23', 'RDFV20', 'RDFV21', 'L2e', 'MoRSEV5', 'RDFV28', 'RDFV29',
       'MoRSEP3', 'P1u', 'rygr', 'Ve', 'MoRSEE7', 'MoRSEV4', 'MoRSEP2', 'FrTATP', 'MoRSEE5', 'P2v', 'ASA', 'MoRSEC21',
       'MoRSEV3', 'MoRSEE1', 'E3v', 'E3u', 'Ke', 'E3p', 'Km', 'MoRSEV7', 'E3e', 'Kp', 'Kv', 'Ku', 'MoRSEV2', 'RDFC8',
       'E3m', 'RDFC9', 'MoRSEM21', 'MoRSEM20', 'MoRSEM23', 'MoRSEM22', 'MoRSEM25', 'MoRSEM24', 'MoRSEM27', 'MoRSEM26',
       'MoRSEM29', 'MoRSEM28', 'MoRSEN19', 'MoRSEN18', 'L2m', 'MoRSEV9', 'MoRSEV8', 'MoRSEU29', 'MoRSEU28', 'L2u',
       'MoRSEV1', 'MoRSEN12', 'RDFC5', 'MoRSEU21', 'MoRSEU20', 'MoRSEU23', 'MoRSEU22', 'MoRSEU25', 'MoRSEU24',
       'MoRSEU27', 'MoRSEU26', 'RDFC7', 'MoRSEV25', 'Tv', 'Am', 'SEig', 'Tu', 'Tp', 'RDFC1', 'Tm', 'RDFC2', 'PNSA3',
       'PNSA2', 'PNSA1', 'RDFC3', 'MSA', 'MoRSEE22', 'L1m', 'De', 'MoRSEE10', 'MoRSEE11', 'MoRSEE12', 'MoRSEE13',
       'MoRSEP9', 'MoRSEE15', 'MoRSEE16', 'MoRSEE17', 'MoRSEE18', 'MoRSEE19', 'Du', 'Dp', 'Vu', 'P2e', 'E1p', 'E1u',
       'E1v', 'E1m', 'RNCS', 'MoRSEP28', 'MoRSEN30', 'E1e', 'MoRSEN7', 'W3D', 'RDFU18', 'RDFU19', 'RDFU14', 'RDFU15',
       'RDFU16', 'RDFU17', 'RDFU10', 'RDFU11', 'RDFU12', 'RDFU13', 'Ap', 'Au', 'RDFC4', 'MoRSEP29', 'Av', 'L3m',
       'RDFM30', 'MoRSEE14', 'E2m']


####################
# COMPUTE PARSING  #
####################



################################################################################
class Atom:
    """
    #################################################################
    A atom class used for wrapping some properties of atoms.

    Note that Coordinates is the output of the function

    (_ReadCoordinates).
    #################################################################
    """

    def __init__(self, Coordinates):

        self.pos = vector3d.Vector3d()
        self.radius = 0.0
        self.Coordinates = Coordinates
        self.Element = ''

    def SetCoordinates(self):

        temp = self.Coordinates
        self.pos.x = float(temp[0])
        self.pos.y = float(temp[1])
        self.pos.z = float(temp[2])

    def GetCoordinates(self):

        self.SetCoordinates()

        return self.pos

    def SetElement(self):

        temp = self.Coordinates

        self.Element = temp[3]

    def GetElement(self):

        self.SetElement()

        return self.Element

    def SetRadius(self):

        radii = {'H': 1.20, 'N': 1.55, 'Na': 2.27, 'Cu': 1.40, 'Cl': 1.75, 'C': 1.70,
                 'O': 1.52, 'I': 1.98, 'P': 1.80, 'B': 1.85, 'Br': 1.85, 'S': 1.80, 'Se': 1.90,
                 'F': 1.47, 'Fe': 1.80, 'K': 2.75, 'Mn': 1.73, 'Mg': 1.73, 'Zn': 1.39, 'Hg': 1.8,
                 'Li': 1.8, '.': 1.8}

        temp = self.GetElement()

        if temp in radii.keys():
            self.radius = radii[temp]
        else:
            self.radius = radii['.']

    def GetRadius(self):

        self.SetRadius()

        return self.radius




###########################################################################

def GetAtomClassList(Coordinates):
    """
    #################################################################
    Combine all atoms in a molecule into a list form.
    Note that Coordinates is the output of the function (_ReadCoordinates).
    #################################################################
    """
    Atoms = []
    for i in Coordinates:
        atom = Atom(i)
        atom.SetCoordinates()
        atom.SetElement()
        atom.SetRadius()
        Atoms.append(atom)
    return Atoms


def GetAtomCoordinateMatrix(lcoordinates):
    """
    #################################################################
    Get the atom coordinate matrix
    #################################################################
    """
    nAtom = len(lcoordinates)
    CoordinateMatrix = scipy.zeros([nAtom, 3])
    AtomLabel = []

    for i, j in enumerate(lcoordinates):
        CoordinateMatrix[i, :] = [j[0], j[1], j[2]]
        AtomLabel.append(j[3])

    return scipy.matrix(CoordinateMatrix), AtomLabel



def parseSDFfor3D(pfilin):
    """
    Read the coordinates and charge of each atom in molecule from .sdf file.
    """

    dchargeSDF = {7:-3.0, 6:-2.0, 5:-1.0, 0:0.0, 3:1.0, 2:2.0, 1:3.0} # and 4 for radical


    latoms = []

    filin = file(pfilin, 'r')
    llines = filin.readlines()
    filin.close()

    # start at line 5 classical format
    for AtBlock in llines[4:]:
        if len(AtBlock) != 70 and len(AtBlock) != 52:
            break
        else:
            #print "-" + AtBlock[0:10] + "-"
            #print AtBlock
            # remove IND from 3D and protonation issues
            if search("IND", AtBlock):
                #print "INNN"
                continue
            X = float(AtBlock[0:10].replace(" ", ""))
            Y = float(AtBlock[10:20].replace(" ", ""))
            Z = float(AtBlock[20:30].replace(" ", ""))
            elem = AtBlock[31:34].replace(" ", "")
            charge = int(AtBlock[36:39].replace(" ", ""))
            charge = dchargeSDF[charge]

            at = [X, Y, Z, elem, charge]
            #print at
            latoms.append(at)

    return latoms


def get_atomicMass(element):

    atomicMass={'H': 1.0079, 'N': 14.0067, 'Na': 22.9897, 'Cu': 63.546, 'Cl': 35.453, 'C': 12.0107,
                 'O': 15.9994, 'I': 126.9045, 'P': 30.9738, 'B': 10.811, 'Br': 79.904, 'S': 32.065, 'Se': 78.96,
                 'F': 18.9984, 'Fe': 55.845, 'K': 39.0983, 'Mn': 54.938, 'Mg': 24.305, 'Zn': 65.39, 'Hg': 200.59,
                 'Li': 6.941, 'Co': 58.9332, "Si":28.0855, "As": 74.9216, "Te":127.6, "Sr":87.62}


    return atomicMass[element]


def get_MW(lcoords, H=1):
    MW = 0
    for coords in lcoords:
        if coords[3] != "H" and H ==0:MW = MW + get_atomicMass(coords[3])
        else:MW = MW + get_atomicMass(coords[3])
    return MW






def get3Ddesc(psdf, geometry=1, cpsa=1, rdf=1, morse=1, whim=1):

    ddesc = {}
    lcoordinates = parseSDFfor3D(psdf)
    print psdf

    if not path.exists(psdf):
        return {}

    if geometry == 1:
        ddesc['W3DH'] = geo3D.Calculate3DWienerWithH(lcoordinates)
        ddesc['W3D'] = geo3D.Calculate3DWienerWithoutH(lcoordinates)
        ddesc['Petitj3D'] = geo3D.CalculatePetitjean3DIndex(lcoordinates)
        ddesc['GeDi'] = geo3D.CalculateGemetricalDiameter(lcoordinates)
        ddesc['grav'] = geo3D.CalculateGravitational3D1(lcoordinates)
        ddesc['rygr'] = geo3D.CalculateRadiusofGyration(lcoordinates)
        ddesc['Harary3D'] = geo3D.CalculateHarary3D(lcoordinates)
        ddesc['AGDD'] = geo3D.CalculateAverageGeometricalDistanceDegree(lcoordinates)
        ddesc['SEig'] = geo3D.CalculateAbsEigenvalueSumOnGeometricMatrix(lcoordinates)


        ddesc['SPAN'] = geo3D.CalculateSPANR(lcoordinates)
        ddesc['ASPAN'] = geo3D.CalculateAverageSPANR(lcoordinates)
        ddesc['MEcc'] = geo3D.CalculateMolecularEccentricity(lcoordinates)

    if cpsa == 1:
        ChargeSA = cpsa3D.GetChargeSA(lcoordinates, RadiusProbe=1.5, n_sphere_point=5000)

        ddesc['ASA'] = cpsa3D.CalculateASA(ChargeSA)
        ddesc['MSA'] = cpsa3D.CalculateMSA(lcoordinates)
        ddesc['PNSA1'] = cpsa3D.CalculatePNSA1(ChargeSA)
        ddesc['PNSA2'] = cpsa3D.CalculatePNSA2(ChargeSA)
        ddesc['PNSA3'] = cpsa3D.CalculatePNSA3(ChargeSA)
        ddesc['PPSA1'] = cpsa3D.CalculatePPSA1(ChargeSA)
        ddesc['PPSA2'] = cpsa3D.CalculatePPSA2(ChargeSA)
        ddesc['PPSA3'] = cpsa3D.CalculatePPSA3(ChargeSA)
        ddesc['DPSA1'] = cpsa3D.CalculateDPSA1(ChargeSA)
        ddesc['DPSA2'] = cpsa3D.CalculateDPSA2(ChargeSA)
        ddesc['DPSA3'] = cpsa3D.CalculateDPSA3(ChargeSA)
        ddesc['FNSA1'] = cpsa3D.CalculateFNSA1(ChargeSA)
        ddesc['FNSA2'] = cpsa3D.CalculateFNSA2(ChargeSA)
        ddesc['FNSA3'] = cpsa3D.CalculateFNSA3(ChargeSA)
        ddesc['FPSA1'] = cpsa3D.CalculateFPSA1(ChargeSA)
        ddesc['FPSA2'] = cpsa3D.CalculateFPSA2(ChargeSA)
        ddesc['FPSA3'] = cpsa3D.CalculateFPSA3(ChargeSA)
        ddesc['WNSA1'] = cpsa3D.CalculateWNSA1(ChargeSA)
        ddesc['WNSA2'] = cpsa3D.CalculateWNSA2(ChargeSA)
        ddesc['WNSA3'] = cpsa3D.CalculateWNSA3(ChargeSA)
        ddesc['WPSA1'] = cpsa3D.CalculateWPSA1(ChargeSA)
        ddesc['WPSA2'] = cpsa3D.CalculateWPSA2(ChargeSA)
        ddesc['WPSA3'] = cpsa3D.CalculateWPSA3(ChargeSA)
        ddesc['TASA'] = cpsa3D.CalculateTASA(ChargeSA)
        ddesc['PSA'] = cpsa3D.CalculateTPSA(ChargeSA)
        ddesc['RASA'] = cpsa3D.CalculateRASA(ChargeSA)
        ddesc['RPSA'] = cpsa3D.CalculateRPSA(ChargeSA)
        ddesc['RNCS'] = cpsa3D.CalculateRNCS(ChargeSA)
        ddesc['RPCS'] = cpsa3D.CalculateRPCS(ChargeSA)
        ddesc['FrTATP'] = cpsa3D.CalculateFractionTATP(ChargeSA)

    if rdf ==1:
        ddesc.update(rdf3D.CalculateUnweightRDF(lcoordinates))
        ddesc.update(rdf3D.CalculateChargeRDF(lcoordinates))
        ddesc.update(rdf3D.CalculateMassRDF(lcoordinates))
        ddesc.update(rdf3D.CalculatePolarizabilityRDF(lcoordinates))
        ddesc.update(rdf3D.CalculateSandersonElectronegativityRDF(lcoordinates))
        ddesc.update(rdf3D.CalculateVDWVolumeRDF(lcoordinates))

    if morse ==1:
        ddesc.update(morse3D.CalculateUnweightMoRSE(lcoordinates))
        ddesc.update(morse3D.CalculateChargeMoRSE(lcoordinates))
        ddesc.update(morse3D.CalculateMassMoRSE(lcoordinates))
        ddesc.update(morse3D.CalculateAtomicNumberMoRSE(lcoordinates))
        ddesc.update(morse3D.CalculatePolarizabilityMoRSE(lcoordinates))
        ddesc.update(morse3D.CalculateSandersonElectronegativityMoRSE(lcoordinates))
        ddesc.update(morse3D.CalculateVDWVolumeMoRSE(lcoordinates))

    if whim ==1:
        CoordinateMatrix, AtomLabel = GetAtomCoordinateMatrix(lcoordinates)
        ddesc['L1u'] = whim3D.GetWHIM1(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['L2u'] = whim3D.GetWHIM2(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['L3u'] = whim3D.GetWHIM3(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['Tu'] = whim3D.GetWHIM4(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['Au'] = whim3D.GetWHIM5(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['Vu'] = whim3D.GetWHIM6(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['P1u'] = whim3D.GetWHIM7(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['P2u'] = whim3D.GetWHIM8(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['Ku'] = whim3D.GetWHIM9(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['E1u'] = whim3D.GetWHIM10(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['E2u'] = whim3D.GetWHIM11(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['E3u'] = whim3D.GetWHIM12(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['Du'] = whim3D.GetWHIM13(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['L1m'] = whim3D.GetWHIM1(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['L2m'] = whim3D.GetWHIM2(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['L3m'] = whim3D.GetWHIM3(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['Tm'] = whim3D.GetWHIM4(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['Am'] = whim3D.GetWHIM5(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['Vm'] = whim3D.GetWHIM6(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['P1m'] = whim3D.GetWHIM7(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['P2m'] = whim3D.GetWHIM8(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['Km'] = whim3D.GetWHIM9(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['E1m'] = whim3D.GetWHIM10(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['E2m'] = whim3D.GetWHIM11(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['E3m'] = whim3D.GetWHIM12(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['Dm'] = whim3D.GetWHIM13(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['L1e'] = whim3D.GetWHIM1(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['L2e'] = whim3D.GetWHIM2(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['L3e'] = whim3D.GetWHIM3(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['Te'] = whim3D.GetWHIM4(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['Ae'] = whim3D.GetWHIM5(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['Ve'] = whim3D.GetWHIM6(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['P1e'] = whim3D.GetWHIM7(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['P2e'] = whim3D.GetWHIM8(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['Ke'] = whim3D.GetWHIM9(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['E1e'] = whim3D.GetWHIM10(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['E2e'] = whim3D.GetWHIM11(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['E3e'] = whim3D.GetWHIM12(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['De'] = whim3D.GetWHIM13(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['L1v'] = whim3D.GetWHIM1(CoordinateMatrix, AtomLabel, proname='V')
        try:ddesc['L2v'] = whim3D.GetWHIM2(CoordinateMatrix, AtomLabel, proname='V')
        except: ddesc['L2v']="NA"
        ddesc['L3v'] = whim3D.GetWHIM3(CoordinateMatrix, AtomLabel, proname='V')
        ddesc['Tv'] = whim3D.GetWHIM4(CoordinateMatrix, AtomLabel, proname='V')
        ddesc['Av'] = whim3D.GetWHIM5(CoordinateMatrix, AtomLabel, proname='V')
        ddesc['Vv'] = whim3D.GetWHIM6(CoordinateMatrix, AtomLabel, proname='V')
        ddesc['P1v'] = whim3D.GetWHIM7(CoordinateMatrix, AtomLabel, proname='V')
        ddesc['P2v'] = whim3D.GetWHIM8(CoordinateMatrix, AtomLabel, proname='V')
        ddesc['Kv'] = whim3D.GetWHIM9(CoordinateMatrix, AtomLabel, proname='V')
        ddesc['E1v'] = whim3D.GetWHIM10(CoordinateMatrix, AtomLabel, proname='V')
        ddesc['E2v'] = whim3D.GetWHIM11(CoordinateMatrix, AtomLabel, proname='V')
        ddesc['E3v'] = whim3D.GetWHIM12(CoordinateMatrix, AtomLabel, proname='V')
        ddesc['Dv'] = whim3D.GetWHIM13(CoordinateMatrix, AtomLabel, proname='V')
        ddesc['L1p'] = whim3D.GetWHIM1(CoordinateMatrix, AtomLabel, proname='alapha')
        ddesc['L2p'] = whim3D.GetWHIM2(CoordinateMatrix, AtomLabel, proname='alapha')
        ddesc['L3p'] = whim3D.GetWHIM3(CoordinateMatrix, AtomLabel, proname='alapha')
        ddesc['Tp'] = whim3D.GetWHIM4(CoordinateMatrix, AtomLabel, proname='alapha')
        ddesc['Ap'] = whim3D.GetWHIM5(CoordinateMatrix, AtomLabel, proname='alapha')
        ddesc['Vp'] = whim3D.GetWHIM6(CoordinateMatrix, AtomLabel, proname='alapha')
        ddesc['P1p'] = whim3D.GetWHIM7(CoordinateMatrix, AtomLabel, proname='alapha')
        ddesc['P2p'] = whim3D.GetWHIM8(CoordinateMatrix, AtomLabel, proname='alapha')
        ddesc['Kp'] = whim3D.GetWHIM9(CoordinateMatrix, AtomLabel, proname='alapha')
        ddesc['E1p'] = whim3D.GetWHIM10(CoordinateMatrix, AtomLabel, proname='alapha')
        ddesc['E2p'] = whim3D.GetWHIM11(CoordinateMatrix, AtomLabel, proname='alapha')
        ddesc['E3p'] = whim3D.GetWHIM12(CoordinateMatrix, AtomLabel, proname='alapha')
        ddesc['Dp'] = whim3D.GetWHIM13(CoordinateMatrix, AtomLabel, proname='alapha')
        ddesc['P3p'] = whim3D.GetWHIM14(CoordinateMatrix, AtomLabel, proname='alapha')
        ddesc['P3u'] = whim3D.GetWHIM14(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['P3m'] = whim3D.GetWHIM14(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['P3e'] = whim3D.GetWHIM14(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['P3v'] = whim3D.GetWHIM14(CoordinateMatrix, AtomLabel, proname='V')

    return ddesc


