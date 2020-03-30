from numpy import *
from copy import deepcopy


def groupAtomCoord(list_atom1):
    l_out = []
    for atom in list_atom1:
        l_out.append([float(atom.x), float(atom.y), float(atom.z)])

    return l_out


# def transformMatrix2List (mat_cord):
#
#     l_out = []
#     group_atom_rotated = array(mat_cord) # tranform array
#     for rowmat in group_atom_rotated :
#         d_atom = {}
#         d_atom["x"] = rowmat[0]
#         d_atom["y"] = rowmat[1]
#         d_atom["z"] = rowmat[2]
#         d_atom["element"] = "C"
#         d_atom["resSeq"] = "LIZ"
#
#         l_out.append(d_atom)
#     return l_out


def rigid_transform_3D(l_at1, l_at2):
    # tranformation array
    A = mat(array(groupAtomCoord(l_at1)))
    B = mat(array(groupAtomCoord(l_at2)))

    if len(A) != len(B):
        print "Not same number of atoms"
        return None, None
    assert len(A) == len(B)

    N = A.shape[0]  # total points

    centroid_A = sum(A, axis=0) / float(N)
    centroid_B = sum(B, axis=0) / float(N)

    #     # centre the points
    AA = A - tile(centroid_A, (N, 1))
    #     print tile(centroid_A, (N, 1)), "11111111"
    BB = B - tile(centroid_B, (N, 1))
    #
    #     # dot is matrix multiplication for array
    H = transpose(AA) * BB
    #
    U, S, Vt = linalg.svd(H)
    #
    R = Vt.T * U.T
    t = -R * centroid_A.T + centroid_B.T
    #
    #     # special reflection case
    if linalg.det(R) < 0:
        print "Reflection detected"
        Vt[2, :] *= -1
        R2 = Vt.T * U.T
        t2 = -R2 * centroid_A.T + centroid_B.T

    if "R2" in locals():
        A1 = applyTranformation(R, t, A)
        A2 = applyTranformation(R2, t2, A)
        rmse1 = rmse(A1, B)
        rmse2 = rmse(B, A2)
        if rmse2 < rmse1:
            return R2, t2


    return R, t


def applyTranformation(matrix_rotation, matrix_transloc, matrix_points=None, l_atom_in=[]):
    # case with list of point in input
    if matrix_points == None:
        matrix_points = mat(array(groupAtomCoord(l_atom_in)))

    # A2 new coord
    n = len(matrix_points)  # number of points
    A2 = (matrix_rotation * matrix_points.T) + tile(matrix_transloc, (1, n))
    A2 = A2.T

    # return only matrix if list atom in input -> empty

    if l_atom_in == []:
        return A2
    else:
        return transformListPoint(l_atom_in, A2)


def transformListPoint(l_atom_in, matrix_coord_before_transloc):
    l_out = deepcopy(l_atom_in)
    nb_atom = len(matrix_coord_before_transloc)

    i = 0
    while i < nb_atom:
        l_out[i].x = matrix_coord_before_transloc[i].tolist()[0][0]
        l_out[i].y = matrix_coord_before_transloc[i].tolist()[0][1]
        l_out[i].z = matrix_coord_before_transloc[i].tolist()[0][2]

        i = i + 1

    return l_out


def rmse(points1, points2):
    if type(points1) == list:
        m_points2 = mat(array(groupAtomCoord(points2)))
        m_points1 = mat(array(groupAtomCoord(points1)))
    else:
        m_points2 = points2
        m_points1 = points1

    # print m_points1
    #     print m_points2

    # tranform array
    # condition same size matrix
    n = len(m_points1)
    err = m_points2 - m_points1
    err = multiply(err, err)
    err = sum(err)
    f_rmse = sqrt(err / n)

    return f_rmse
