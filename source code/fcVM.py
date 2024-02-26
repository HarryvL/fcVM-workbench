# ***************************************************************************
# *                                                                         *
# *   Copyright (c) 2024 - Harry van Langen <hvlanalysis@gmail.com>        *
# *                                                                         *
# *   This program is free software; you can redistribute it and/or modify  *
# *   it under the terms of the GNU Lesser General Public License (LGPL)    *
# *   as published by the Free Software Foundation; either version 2 of     *
# *   the License, or (at your option) any later version.                   *
# *   for detail see the LICENCE text file.                                 *
# *                                                                         *
# *   This program is distributed in the hope that it will be useful,       *
# *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
# *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
# *   GNU Library General Public License for more details.                  *
# *                                                                         *
# *   You should have received a copy of the GNU Library General Public     *
# *   License along with this program; if not, write to the Free Software   *
# *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  *
# *   USA                                                                   *
# *                                                                         *
# ***************************************************************************


import time
import FemGui
import ObjectsFem
import numpy as np
import Part as Part
import FreeCAD as App
import FreeCADGui as Gui
from FreeCAD import Units
import scipy.sparse as scsp
from numba import jit, types
from numba.typed import Dict
import matplotlib.pyplot as plt
from femtools import membertools
from femmesh import meshsetsgetter
from femmesh import meshtools as mt
from sksparse.cholmod import cholesky
from matplotlib.widgets import Button
from femresult import resulttools as rt
from feminout import importToolsFem as itf

np.set_printoptions(precision=5, linewidth=300)


def prn_upd(*args):
    for obj in args:
        print(str(obj), end='')
    print('\n')
    Gui.updateGui()


def setUpAnalysis():
    doc = App.ActiveDocument
    mesh = doc.getObject("FEMMeshGmsh").FemMesh
    if mesh is None:
        prn_upd("No Gmsh object. Please create one first")
        raise SystemExit()
    analysis = doc.getObject("Analysis")
    if analysis is None:
        prn_upd("No Analysis object. Please create one first")
        raise SystemExit()
    # purge result objects
    rt.purge_results(analysis)
    doc.recompute()

    return doc, mesh, analysis


def setUpInput(doc, mesh, analysis):
    analysis = doc.getObject("Analysis")
    solver = doc.getObject("SolverCcxTools")
    docmesh = doc.getObject("FEMMeshGmsh")
    member = membertools.AnalysisMember(analysis)

    if solver == None:
        FemGui.setActiveAnalysis(App.activeDocument().Analysis)
        FemGui.getActiveAnalysis().addObject(ObjectsFem.makeSolverCalculixCcxTools(App.ActiveDocument))
        solver = doc.getObject("SolverCcxTools")

    # determine elements connected to a node using FC API
    fet = mt.get_femelement_table(mesh)
    # fet is dictionary: { elementid : [ nodeid, nodeid, ... , nodeid ] }
    net = mt.get_femnodes_ele_table(mesh.Nodes, fet)
    # net is dictionary: {nodeID : [[eleID, binary node position], [], ...], nodeID : [[], [], ...], ...}
    # node0 has binary node position 2^0 = 1, node1 = 2^1 = 2, ..., node10 = 2^10 = 1024

    # create connectivity array elNodes for mapping local node number -> global node number
    elNodes = np.array([mesh.getElementNodes(el) for el in mesh.Volumes])  # elNodes[elementIndex] = [node1,...,Node10]

    # create nodal coordinate array nocoord for node number -> (x,y,z)
    ncv = list(mesh.Nodes.values())
    nocoord = np.asarray([[v.x, v.y, v.z] for v in ncv])  # nocoord[nodeIndex] = [x-coord, y-coord, z-coord]

    # get access to element sets: meshdatagetter.mat_geo_sets
    meshdatagetter = meshsetsgetter.MeshSetsGetter(
        analysis,
        solver,
        docmesh,
        member)
    meshdatagetter.get_mesh_sets()

    if len(member.mats_linear) == 1:
        element_sets = [mesh.Volumes]
    else:
        element_sets = [es["FEMElements"] for es in member.mats_linear]

    matCon = {}  # BooleanFragment Primitive the material object refers to
    ppEl = {}  # BooleanFragment Primitive element El belongs to
    materialbyElement = []  # see further

    prn_upd("Number of material objects: ", len(member.mats_linear))

    for indm, matobject in enumerate(member.mats_linear):
        E = float(App.Units.Quantity(matobject['Object'].Material['YoungsModulus']).getValueAs('MPa'))
        Nu = float(matobject['Object'].Material['PoissonRatio'])
        Density = float(App.Units.Quantity(matobject['Object'].Material['Density']).getValueAs('kg/mm^3'))
        prn_upd("Material Object: ", matobject['Object'].Name, "   E= ", E, "   Nu= ", Nu, "   Density= ", Density)
        for el in element_sets[indm]:  # element_sets[indm]: all elements with material indm
            if matCon: ppEl[el] = matCon[indm]  # ppEl[el]: primitive el belongs to
            materialbyElement.append([E, Nu, Density])  # materialbyElement[elementIndex] = [E, Nu, Density]

    materialbyElement = np.asarray(materialbyElement)

    noce = np.zeros((len(nocoord)), dtype=np.int16)
    for i in net:
        noce[i - 1] = len(net[i])

    # create boundary condition array dispfaces
    u0 = Units.Quantity(0.0, 1)
    dispfaces = []
    for obj in App.ActiveDocument.Objects:
        if obj.isDerivedFrom('Fem::ConstraintFixed') or obj.isDerivedFrom('Fem::ConstraintDisplacement'):

            bcnodes = []

            if obj.isDerivedFrom('Fem::ConstraintFixed'):
                bctype = [False, False, False]
                bcvalue = [u0, u0, u0]
            else:
                bctype = [obj.xFree, obj.yFree, obj.zFree]
                bcvalue = [obj.xDisplacement, obj.yDisplacement, obj.zDisplacement]

            for part, boundaries in obj.References:
                for boundary in boundaries:
                    ref = part.Shape.getElement(boundary)
                    if type(ref) == Part.Vertex:
                        bc = mesh.getNodesByVertex(ref)
                        for bcn in bc: bcnodes.append(bcn)
                    elif type(ref) == Part.Edge:
                        bc = mesh.getNodesByEdge(ref)
                        for bcn in bc: bcnodes.append(bcn)
                    elif type(ref) == Part.Face:
                        bc = mesh.getNodesByFace(
                            ref)  # all nodes on a primitive face with a displacement boundary condition
                        for bcn in bc:
                            bcnodes.append(bcn)
                    else:
                        prn_upd("No Boundaries Found")
            bcnodes = list(dict.fromkeys(bcnodes))  # remove duplicates in bcnodes
            if bcnodes:
                dispfaces.append([bcnodes, bctype, bcvalue])

    fix = Dict.empty(key_type=types.int64, value_type=types.float64)
    nn = 3 * len(nocoord)
    fixdof = np.ones(nn, dtype=int)
    movdof = np.zeros(nn, dtype=int)

    for face in dispfaces:
        if not face[1][0]:
            for node in face[0]:
                dof = 3 * (node - 1)
                fix[dof] = face[2][0].Value
                fixdof[dof] = 0
        if not face[1][1]:
            for node in face[0]:
                dof = 3 * (node - 1) + 1
                fix[dof] = face[2][1].Value
                fixdof[dof] = 0
        if not face[1][2]:
            for node in face[0]:
                dof = 3 * (node - 1) + 2
                fix[dof] = face[2][2].Value
                fixdof[dof] = 0

    for dof in fix:
        if fix[dof] != 0.0:
            movdof[dof] = 1

    lf = [[0, 0, 0, 0, 0, 0]]  # load face nodes - signature for numba
    pr = [0.0]  # load face pressure - signature for numba

    for obj in App.ActiveDocument.Objects:
        if obj.isDerivedFrom('Fem::ConstraintPressure'):
            if obj.Reversed:
                sign = 1
            else:
                sign = -1
            for part, faces in obj.References:  # obj.References: references to loaded primitive faces
                for face in faces:
                    ref = part.Shape.getElement(face)
                    if type(ref) == Part.Face:
                        for faceID in mesh.getFacesByFace(ref):  # face ID: ID of a 6-node face element
                            face_nodes = list(mesh.getElementNodes(faceID))  # 6-node element node numbers
                            lf.append(face_nodes)
                            pr.append(sign * obj.Pressure)
                    else:
                        prn_upd("No Faces with Pressure Loads")

    loadfaces = np.array(lf)
    pressure = np.array(pr)

    # re-order element nodes
    for el in elNodes:
        el[1], el[2] = el[2], el[1]
        el[4], el[6] = el[6], el[4]
        el[8], el[9] = el[9], el[8]

    return (elNodes, nocoord, fix, fixdof, movdof, loadfaces, materialbyElement, noce, pressure)


# shape functions for a 4-node tetrahedron - only used for stress interpolation
def shape4tet(xi, et, ze, xl):
    shp = np.zeros((4), dtype=np.float64)

    # shape functions
    shp[0] = 1.0 - xi - et - ze
    shp[1] = xi
    shp[2] = et
    shp[3] = ze

    return shp


@jit(nopython=True, cache=True)
def shp10tet(xi, et, ze):
    shp = np.zeros((10), dtype=np.float64)

    # shape functions - source: Calculix, G Dhondt
    a = 1.0 - xi - et - ze
    shp[0] = (2.0 * a - 1.0) * a
    shp[1] = xi * (2.0 * xi - 1.0)
    shp[2] = et * (2.0 * et - 1.0)
    shp[3] = ze * (2.0 * ze - 1.0)
    shp[4] = 4.0 * xi * a
    shp[5] = 4.0 * xi * et
    shp[6] = 4.0 * et * a
    shp[7] = 4.0 * ze * a
    shp[8] = 4.0 * xi * ze
    shp[9] = 4.0 * et * ze
    return shp


@jit("float64(float64, float64, float64, float64[:,:], float64[:,:])", nopython=True, cache=True, fastmath=True)
def dshp10tet(xi, et, ze, xl, bmat):
    dshp = np.zeros((3, 10), dtype=np.float64)
    dshpg = np.zeros((3, 10), dtype=np.float64)
    xs = np.zeros((3, 3), dtype=np.float64)
    xsi = np.zeros((3, 3), dtype=np.float64)

    # local derivatives of the shape functions: xi-derivative - source: Calculix, G Dhondt
    dshp[0][0] = 1.0 - 4.0 * (1.0 - xi - et - ze)
    dshp[0][1] = 4.0 * xi - 1.0
    # dshp[0][2] = 0.0
    # dshp[0][3] = 0.0
    dshp[0][4] = 4.0 * (1.0 - 2.0 * xi - et - ze)
    dshp[0][5] = 4.0 * et
    dshp[0][6] = -4.0 * et
    dshp[0][7] = -4.0 * ze
    dshp[0][8] = 4.0 * ze
    # dshp[0][9] = 0.0

    # local derivatives of the shape functions: eta-derivative - source: Calculix, G Dhondt
    dshp[1][0] = 1.0 - 4.0 * (1.0 - xi - et - ze)
    # dshp[1][1] = 0.0
    dshp[1][2] = 4.0 * et - 1.0
    # dshp[1][3] = 0.0
    dshp[1][4] = -4.0 * xi
    dshp[1][5] = 4.0 * xi
    dshp[1][6] = 4.0 * (1.0 - xi - 2.0 * et - ze)
    dshp[1][7] = -4.0 * ze
    # dshp[1][8] = 0.0
    dshp[1][9] = 4.0 * ze

    # local derivatives of the shape functions: zeta-derivative - source: Calculix, G Dhondt
    dshp[2][0] = 1.0 - 4.0 * (1.0 - xi - et - ze)
    # dshp[2][1] = 0.0
    # dshp[2][2] = 0.0
    dshp[2][3] = 4.0 * ze - 1.0
    dshp[2][4] = -4.0 * xi
    # dshp[2][5] = 0.0
    dshp[2][6] = -4.0 * et
    dshp[2][7] = 4.0 * (1.0 - xi - et - 2.0 * ze)
    dshp[2][8] = 4.0 * xi
    dshp[2][9] = 4.0 * et

    # xs = np.dot(xl, dshp.T) # local derivative of the global coordinates

    for i in range(3):
        for j in range(3):
            xs[i][j] = 0.0
            for k in range(10):
                xs[i][j] += xl[k][i] * dshp[j][k]

    # xsj = np.linalg.det(xs) # Jacobian

    xsj = (xs[0][0] * xs[1][1] * xs[2][2] -
           xs[0][0] * xs[1][2] * xs[2][1] +
           xs[0][2] * xs[1][0] * xs[2][1] -
           xs[0][2] * xs[1][1] * xs[2][0] +
           xs[0][1] * xs[1][2] * xs[2][0] -
           xs[0][1] * xs[1][0] * xs[2][2])

    # xsi = np.linalg.inv(xs) # global derivative of the local coordinates

    xsi[0][0] = (xs[1][1] * xs[2][2] - xs[2][1] * xs[1][2]) / xsj
    xsi[0][1] = (xs[0][2] * xs[2][1] - xs[0][1] * xs[2][2]) / xsj
    xsi[0][2] = (xs[0][1] * xs[1][2] - xs[0][2] * xs[1][1]) / xsj
    xsi[1][0] = (xs[1][2] * xs[2][0] - xs[1][0] * xs[2][2]) / xsj
    xsi[1][1] = (xs[0][0] * xs[2][2] - xs[0][2] * xs[2][0]) / xsj
    xsi[1][2] = (xs[1][0] * xs[0][2] - xs[0][0] * xs[1][2]) / xsj
    xsi[2][0] = (xs[1][0] * xs[2][1] - xs[2][0] * xs[1][1]) / xsj
    xsi[2][1] = (xs[2][0] * xs[0][1] - xs[0][0] * xs[2][1]) / xsj
    xsi[2][2] = (xs[0][0] * xs[1][1] - xs[1][0] * xs[0][1]) / xsj

    # dshp = np.dot(xsi.T, dshp) # global derivatives of the shape functions

    for i in range(3):
        for j in range(10):
            for k in range(3):
                dshpg[i][j] += xsi[k][i] * dshp[k][j]

    # computation of the strain interpolation matrix bmat
    for i in range(10):
        i3 = 3 * i
        d00 = dshpg[0][i]
        d10 = dshpg[1][i]
        d20 = dshpg[2][i]
        bmat[0][i3] = d00
        bmat[1][i3 + 1] = d10
        bmat[2][i3 + 2] = d20
        bmat[3][i3] = d10
        bmat[3][i3 + 1] = d00
        bmat[4][i3] = d20
        bmat[4][i3 + 2] = d00
        bmat[5][i3 + 1] = d20
        bmat[5][i3 + 2] = d10

    return xsj


# shape functions and their derivatives for a 6-node triangular interface element
@jit(nopython=True, cache=True)
def shape6tri(xi, et, xl):
    shp = np.zeros((6), dtype=np.float64)
    dshp = np.zeros((2, 6), dtype=np.float64)
    bmat = np.zeros((3, 36), dtype=np.float64)

    # shape functions
    shp[0] = (1.0 - xi - et) * (1.0 - 2.0 * xi - 2.0 * et)
    shp[1] = xi * (2.0 * xi - 1.0)
    shp[2] = et * (2.0 * et - 1.0)
    shp[3] = 4.0 * xi * (1.0 - xi - et)
    shp[4] = 4.0 * xi * et
    shp[5] = 4.0 * et * (1 - xi - et)

    # local derivatives of the shape functions: xi-derivative
    dshp[0][0] = -3.0 + 4.0 * et + 4.0 * xi
    dshp[0][1] = -1.0 + 4.0 * xi
    dshp[0][2] = 0.0
    dshp[0][3] = -4.0 * (-1.0 + et + 2.0 * xi)
    dshp[0][4] = 4.0 * et
    dshp[0][5] = -4.0 * et

    # local derivatives of the shape functions: eta-derivative
    dshp[1][0] = -3.0 + 4.0 * et + 4.0 * xi
    dshp[1][1] = 0.0
    dshp[1][2] = -1.0 + 4.0 * et
    dshp[1][3] = -4.0 * xi
    dshp[1][4] = 4.0 * xi
    dshp[1][5] = -4.0 * (-1.0 + 2.0 * et + xi)

    xs = np.dot(dshp, xl.T)  # xs = [ [[dx/dxi],[dy/dxi],[dz/dxi]] , [[dx/det],[dy/det],[dz/det]] ]

    xp = np.cross(xs[0], xs[1])  # vector normal to surface

    xsj = np.linalg.norm(xp)  # Jacobian

    # xsj = np.sqrt(xp[0]*xp[0]+xp[1]*xp[1]+xp[2]*xp[2])

    xx = xs[0] / np.linalg.norm(xs[0])  # unit vector in xi direction

    # xx = np.sqrt(xs[0]*xs[0]+xs[1]*xs[1]+xs[2]*xs[2])

    xp /= xsj  # unit vector normal to surface
    xt = np.cross(xp, xx)  # unit vector tangential to surface and normal to xx

    # computation of the "strain" interpolation matrix bmat
    for i in range(6):
        ia = 3 * i
        ib = ia + 18
        ni = shp[i]
        bmat[0][ia] = ni
        bmat[1][ia + 1] = ni
        bmat[2][ia + 2] = ni
        bmat[0][ib] = -ni
        bmat[1][ib + 1] = -ni
        bmat[2][ib + 2] = -ni

    return xsj, shp, bmat, xx, xt, xp


# linear-elastic material stiffness matrix

@jit(nopython=True, cache=True)
def hooke(element, materialbyElement, dmat):
    e = materialbyElement[element][0]  # Young's Modulus
    nu = materialbyElement[element][1]  # Poisson's Ratio
    dm = e * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu)
    od = nu / (1.0 - nu)
    sd = 0.5 * (1.0 - 2.0 * nu) / (1.0 - nu)

    dmat[0][0] = dmat[1][1] = dmat[2][2] = 1.0
    dmat[3][3] = dmat[4][4] = dmat[5][5] = sd
    dmat[0][1] = dmat[0][2] = dmat[1][2] = od
    dmat[1][0] = dmat[2][0] = dmat[2][1] = od
    dmat *= dm


# Gaussian integration points and weights
@jit(nopython=True, cache=True)
def gaussPoints():
    # Gaussian integration points and weights for 10-noded tetrahedron
    gp10 = np.array([[0.138196601125011, 0.138196601125011, 0.138196601125011,
                      0.041666666666667],
                     [0.585410196624968, 0.138196601125011, 0.138196601125011,
                      0.041666666666667],
                     [0.138196601125011, 0.585410196624968, 0.138196601125011,
                      0.041666666666667],
                     [0.138196601125011, 0.138196601125011, 0.585410196624968,
                      0.041666666666667]])
    # Gaussian integration points and weights for 6-noded triangle
    gp6 = np.array([[0.445948490915965, 0.445948490915965,
                     0.111690794839005],
                    [0.10810301816807, 0.445948490915965,
                     0.111690794839005],
                    [0.445948490915965, 0.10810301816807,
                     0.111690794839005],
                    [0.091576213509771, 0.091576213509771,
                     0.054975871827661],
                    [0.816847572980458, 0.091576213509771,
                     0.054975871827661],
                    [0.091576213509771, 0.816847572980458,
                     0.054975871827661]])
    return gp10, gp6


# calculate the global stiffness matrix and load vector
# @jit(
#     "types.Tuple((float64[:],int64[:], int64[:], float64[:], float64[:]))(int64[:,:], float64[:,:], float64[:,:], DictType(int64,float64), int64[:,:], float64, float64, float64, float64[:])",
#     nopython=True, cache=True)
@jit(nopython=True, cache=True)
def calcGSM(elNodes, nocoord, materialbyElement, fix, loadfaces, grav_x, grav_y, grav_z, pressure):
    gp10, gp6 = gaussPoints()
    ne = len(elNodes)  # number of volume elements
    nn = len(nocoord[:, 0])  # number of degrees of freedom

    # number of entries in the lower diagonal of the stiffness matrix: (dof**2 + dof) / 2
    ns = int((30 * 30 + 30) / 2 * ne)

    print(f"ne = {ne}")
    print(f"ns = {ns}\n")

    xlf = np.zeros((3, 6), dtype=types.float64)  # coordinates of load face nodes
    xlv = np.zeros((10, 3), dtype=types.float64)  # coordinates of volume element nodes
    glv = np.zeros((3 * nn), dtype=types.float64)  # global load vector
    row = np.zeros(ns, dtype=types.int64)  # row indices of COO matrix
    col = np.zeros(ns, dtype=types.int64)  # column indices of COO matrix
    stm = np.zeros(ns, dtype=types.float64)  # stiffness values of COO matrix
    modf = np.zeros((3 * nn), dtype=types.float64)  # modification to the stiffness matrix for displacement BCs
    dof = np.zeros(30, dtype=types.int64)
    dmat = np.zeros((6, 6), dtype=np.float64)
    bmatV = np.zeros((6, 30), dtype=np.float64)

    #   calculate element load vectors for pressure and add to global vector

    for face in range(len(pressure) - 1):  # first pressure value is a dummy signature for numba
        if len(pressure) == 1:
            break

        nda = loadfaces[face + 1]  # node numbers of loaded face
        for i in range(3):
            for j in range(6):
                nd = nda[j]
                xlf[i][j] = nocoord[nd - 1][i]  # coordinates of loaded face nodes

        # integrate element load vector
        for index in range(len(gp6)):
            xi = gp6[index][0]
            et = gp6[index][1]
            xsj, shp, bmatS, xx, xt, xp = shape6tri(xi, et, xlf)
            nl = 0
            for i in range(len(loadfaces[face] - 1)):
                nd = loadfaces[face + 1][i]
                iglob = nd - 1
                iglob3 = 3 * iglob
                for k in range(3):
                    load = shp[nl] * pressure[face + 1] * xp[k] * abs(xsj) * gp6[index][2]
                    glv[iglob3 + k] += load
                nl += 1

    # for each volume element calculate the element stiffness matrix
    # and gravity load vector and add to global matrix and vector

    pos = 0

    V = 0.0

    hooke(0, materialbyElement, dmat)
    density = materialbyElement[0][2]

    for el, nodes in enumerate(elNodes):
        esm = np.zeros((30, 30), dtype=np.float64)
        gamma = np.zeros((30), dtype=np.float64)

        # set up nodal values for this element
        for i in range(3):
            for j in range(10):
                nd = nodes[j]
                xlv[j][i] = nocoord[nd - 1][i]

        # integrate element matrix
        for ip in gp10:
            xi = ip[0]
            et = ip[1]
            ze = ip[2]
            shp = shp10tet(xi, et, ze)
            xsj = dshp10tet(xi, et, ze, xlv, bmatV)
            esm += np.dot(bmatV.T, np.dot(dmat, bmatV)) * ip[3] * abs(xsj)
            gamma[0::3] += grav_x * density * shp * ip[3] * abs(xsj)
            gamma[1::3] += grav_y * density * shp * ip[3] * abs(xsj)
            gamma[2::3] += grav_z * density * shp * ip[3] * abs(xsj)
            V += xsj * ip[3]  # Element volume - not used

        for i in range(10):
            nd = nodes[i] - 1
            glv[3 * nd] += gamma[3 * i]
            glv[3 * nd + 1] += gamma[3 * i + 1]
            glv[3 * nd + 2] += gamma[3 * i + 2]
            for j in range(3):
                dof[3 * i + j] = 3 * nd + j

        for i in range(30):
            dofi = dof[i]
            if dofi in fix:
                row[pos] = dofi
                col[pos] = dofi
                stm[pos] = 1.0
                modf[dofi] += fix[dofi]
                pos += 1
                for j in range(i):
                    dofj = dof[j]
                    if dofj not in fix:
                        modf[dofj] -= esm[i][j] * fix[dofi]
            else:
                for j in range(i + 1):
                    dofj = dof[j]
                    if dofj in fix:
                        modf[dofi] -= esm[i][j] * fix[dofj]
                    else:
                        if dofi > dofj:
                            row[pos] = dofi
                            col[pos] = dofj
                        else:
                            row[pos] = dofj
                            col[pos] = dofi
                        stm[pos] = esm[i][j]
                        pos += 1

    row = row[:pos]
    col = col[:pos]
    stm = stm[:pos]

    loadsumx = 0.0
    loadsumy = 0.0
    loadsumz = 0.0
    for node in range(nn):
        dof = 3 * node
        loadsumx += glv[dof]
        loadsumy += glv[dof + 1]
        loadsumz += glv[dof + 2]

    print("mesh volume", V)
    print("loadsumx", loadsumx)
    print("loadsumy", loadsumy)
    print("loadsumz", loadsumz, "\n")

    return stm, row, col, glv, modf


# calculate load-deflection curve
def calcDisp(elNodes, nocoord, fixdof, movdof, modf, materialbyElement, stm, row, col,
             glv, nstep, iterat_max, error_max, relax, scale_re, scale_up, scale_dn, sig_yield, disp_output,
             ultimate_strain, progress_update):
    ndof = len(glv)  # number of degrees of freedom
    nelem = len(elNodes)  # number of elements
    t0 = time.perf_counter()
    gsm = scsp.csc_matrix((stm, (row, col)), shape=(ndof, ndof))  # construct sparse global stiffness matrix
    t1 = time.perf_counter()
    prn_upd("construct sparse global stiffness matrix: {:.2e} s".format((t1 - t0)))

    # TODO: improve error estimate for pure displacement control
    qnorm = np.linalg.norm(glv)
    if qnorm < 1.0: qnorm = 1.0

    # Cholesky decomposition of the global stiffness matrix and elastic solution using Cholmod
    t0 = time.perf_counter()
    factor = cholesky(gsm)
    t1 = time.perf_counter()
    ue = factor(fixdof * glv + modf)  # elastic solution
    t2 = time.perf_counter()
    prn_upd("sparse Cholesky decomposition: {:.2e} s, elastic solution: {:.2e} s".format((t1 - t0), (t2 - t1)))

    # initiate analysis
    dl0 = 1.0 / nstep  # nstep == 1 execute an elastic analysis
    dl = dl0
    du = dl * ue

    sig_new = np.zeros(24 * nelem, dtype=np.float64)  # stress in Tet10
    sig_old = np.zeros(24 * nelem, dtype=np.float64)
    sig_test = np.zeros(24 * nelem, dtype=np.float64)
    peeq = np.zeros(4 * nelem, dtype=np.float64)  # equivalent plastic strain in Tet10
    csr = np.zeros(4 * nelem, dtype=np.float64)  # critical strain ratio in Tet10
    disp_new = np.zeros(ndof, dtype=np.float64)  # displacement results
    disp_old = np.zeros(ndof, dtype=np.float64)  # displacement results
    lbd = np.zeros(1, dtype=np.float64)  # load level
    rfl = np.zeros(1, dtype=np.float64)  # reaction force level (for displacement control)

    gp10, gp6 = gaussPoints()

    # determine elastic reaction force on moving boundary
    if max(movdof) == 1:
        qelastic = np.zeros(3 * len(nocoord), dtype=np.float64)
        update_stress_load(gp10, elNodes, nocoord, materialbyElement, sig_yield, ue, sig_old, sig_new, sig_test,
                           qelastic)

        qelastic *= movdof
        qelnorm = np.linalg.norm(qelastic)
        qnorm = qelnorm
        sig_new = np.zeros(24 * nelem, dtype=np.float64)  # reset sig_new to zero

    step = -1
    cnt = True

    un = [0.]
    csrplot = [0.]

    if float(nstep) == 1.0:
        # perform an elastic (one-step) analysis
        step = 0
        out_disp = 1
        lbd = np.append(lbd, 1.0)
        rfl = np.append(rfl, 1.0)
        disp_new = ue
        un.append(np.max(np.abs(disp_new)))
        cnt = False

    factor_time_tot = 0.0
    iterat_tot = 0

    while cnt:
        cnt = False
        pstep = 0
        progress_update(0)
        for _ in (range(nstep)):
            step += 1
            pstep += 1
            progress_update(int(100 * (pstep + 1) / nstep))
            restart = 0
            sig_old = sig_new.copy()

            prn_upd("Step: {}".format(step))
            a = du.copy()  # a: Riks control vector
            lbd = np.append(lbd, lbd[step] + dl)  # lbd: load level

            # update stresses and loads
            qin = np.zeros(3 * len(nocoord), dtype=np.float64)
            update_stress_load(gp10, elNodes, nocoord, materialbyElement, sig_yield, du, sig_old, sig_new, sig_test,
                               qin)

            # calculate residual load vector
            fex = fixdof * lbd[step + 1] * glv
            fin = fixdof * qin
            r = fex - fin
            rnorm = np.linalg.norm(r)

            # out-of-balance error
            error = rnorm / qnorm

            iterat = 0
            prn_upd("Iteration: {}, Error: {:.2e}".format(iterat, error))

            while error > error_max:

                iterat += 1
                iterat_tot += 1

                t0 = time.perf_counter()
                due = factor(relax * r)
                t1 = time.perf_counter()
                factor_time_tot += t1 - t0

                # Riks control correction to load level increment
                dl = -np.dot(a, due) / np.dot(a, ue)
                lbd[step + 1] += dl

                # Riks control correction to displacement increment
                du += due + dl * ue

                # update stresses and loads

                qin = np.zeros(3 * len(nocoord), dtype=np.float64)
                update_stress_load(gp10, elNodes, nocoord, materialbyElement, sig_yield, du, sig_old, sig_new, sig_test,
                                   qin)

                # calculate out of balance error
                r = fixdof * (lbd[step + 1] * glv - qin)
                rnorm = np.linalg.norm(r)
                error = rnorm / qnorm
                prn_upd("Iteration: {}, Error: {:.2e}".format(iterat, error))

                if iterat > iterat_max:
                    # scale down
                    if restart == 4:
                        print("MAXIMUM RESTARTS REACHED")
                        raise SystemExit()
                    restart += 1
                    if step > 0:
                        dl = (lbd[step] - lbd[step - 1]) / scale_re / restart
                        du = (disp_new - disp_old) / scale_re / restart
                    else:
                        # for first step only
                        dl = dl0 / scale_re / restart
                        du = dl * ue / scale_re / restart
                    lbd[step + 1] = lbd[step] + dl

                    qin = np.zeros(3 * len(nocoord), dtype=np.float64)
                    update_stress_load(gp10, elNodes, nocoord, materialbyElement, sig_yield, du, sig_old, sig_new,
                                       sig_test, qin)

                    r = fixdof * (lbd[step + 1] * (glv + modf) - qin)
                    rnorm = np.linalg.norm(r)
                    error = rnorm / qnorm

                    iterat = 0

            # update results at end of converged load step
            disp_old = disp_new.copy()
            disp_new += du
            dl = lbd[step + 1] - lbd[step]
            if max(movdof) == 1:
                rfl = np.append(rfl, np.sum(movdof * qin))

            if iterat > 10:
                # scale down
                dl /= scale_dn
                du /= scale_dn
            if iterat < 5:
                # scale up
                dl *= scale_up
                du *= scale_up
            un.append(np.max(np.abs(disp_new)))
            update_PEEQ_CSR(nelem, materialbyElement, sig_test, sig_new, sig_yield, ultimate_strain, peeq, csr)
            csrplot.append(np.max(csr))

        csr_non_zero = np.nonzero(csrplot)
        if len(csr_non_zero[0]) != 0:
            el_limit = csr_non_zero[0][0] - 1
        else:
            el_limit = 0
        csr_limit = np.argwhere(np.asarray(csrplot) > 1.0)
        if len(csr_limit) != 0:
            ul_limit = csr_limit[0][0] - 1
        else:
            ul_limit = 0

        # plot load-displacement curve - TODO: move to output / post-processing
        if max(movdof) == 1:
            cnt = plot(el_limit, ul_limit, un, rfl, csrplot)
        else:
            cnt = plot(el_limit, ul_limit, un, lbd, csrplot)

    prn_upd("total time evaluating K_inv * r: {}".format(factor_time_tot))
    prn_upd("total number of iterations: {}".format(iterat_tot))
    if iterat_tot != 0:
        prn_upd("average time to evaluate K_inv * r: {}".format(factor_time_tot / iterat_tot))

    out = step + 1
    u_out = un[out] - un[out - 1]
    l_out = lbd[out] - lbd[out - 1]
    prn_upd("Step: {0:2d} Load level increment: {1:.3f} Displacement increment: {2:.4e}".format(out, l_out,
                                                                                                u_out))
    if disp_output == "total":
        return disp_new, sig_new, peeq
    else:
        return disp_new - disp_old, sig_new, peeq


# plot the load-deflection curve
def plot(el_limit, ul_limit, un, lbd, csrplot):
    class Index(object):
        def stop(self, event):
            self.cnt = False
            plt.close()
            self.clicked = True

        def add(self, event):
            self.cnt = True
            plt.close()
            self.clicked = True

        def close_window(self, event):
            if self.cnt == False:
                self.stop('stop_event')

    callback = Index()
    callback.cnt = False
    callback.clicked = False
    fig, ax = plt.subplots(1, 2, figsize=(10, 6))
    fig.canvas.manager.set_window_title('fcVM')
    plt.subplots_adjust(bottom=0.2)
    ax[0].plot(un, lbd, '-ok')
    ax[0].set(xlabel='displacement [mm]', ylabel='load factor [-] or load [N]', title='')
    ax[1].plot(csrplot, lbd, '-ok')
    ax[1].set(xlabel='critical strain ratio [-]', title='')
    ax[0].grid()
    ax[1].grid()
    axstop = plt.axes([0.5 - 0.075 - 0.0075, 0.05, 0.075, 0.06])
    axadd = plt.axes([0.5 + 0.0075, 0.05, 0.075, 0.06])
    bstop = Button(axstop, 'stop')
    bstop.on_clicked(callback.stop)
    badd = Button(axadd, 'add')
    badd.on_clicked(callback.add)
    fig.canvas.mpl_connect('close_event', callback.close_window)

    if ul_limit != 0:
        fac = (1.0 - csrplot[ul_limit]) / (csrplot[ul_limit + 1] - csrplot[ul_limit])
        lbd_limit = lbd[ul_limit] + fac * (lbd[ul_limit + 1] - lbd[ul_limit])
        un_limit = un[ul_limit] + fac * (un[ul_limit + 1] - un[ul_limit])
        ax[0].plot([0.0, un_limit], [lbd_limit, lbd_limit], color='r', linestyle="--")
        ax[0].plot([un_limit, un_limit], [0.0, lbd_limit], color='r', linestyle="--")
        ax[0].plot([un[el_limit], un[el_limit]], [0.0, lbd[el_limit]], color='b', linestyle="--")
        ax[0].plot([0.0, un[el_limit]], [lbd[el_limit], lbd[el_limit]], color='b', linestyle="--")
        ax[1].plot([0.0, 1.0], [lbd_limit, lbd_limit], color='r', linestyle="--")
        ax[1].plot([1.0, 1.0], [0.0, lbd_limit], color='r', linestyle="--")
        ax[1].plot([0.0, 1.0], [lbd[el_limit], lbd[el_limit]], color='b', linestyle="--")

    plt.show()

    while True:
        plt.pause(0.01)
        if callback.clicked:
            break

    return callback.cnt


# update PEEQ and CSR
@jit(nopython=True, cache=True, fastmath=True)
def update_PEEQ_CSR(nelem, materialbyElement, sig_test, sig_new, sig_yield, ultimate_strain, peeq, csr):
    E = materialbyElement[0][0]  # Young's Modulus
    nu = materialbyElement[0][1]  # Poisson's Ratio
    G = E / 2.0 / (1 + nu)  # shear modulus
    if ultimate_strain == 0.0:
        ultimate_strain = 1.0e12
    alpha = np.sqrt(np.e) * ultimate_strain  # stress triaxiality T = 1/3 for uniaxial test
    beta = 1.5

    for el in range(nelem):
        for ip in range(4):
            ipos1 = 4 * el + ip
            ipos2 = 24 * el + 6 * ip
            st0, st1, st2, st3, st4, st5 = sig_test[ipos2:ipos2 + 6]
            sn0, sn1, sn2, sn3, sn4, sn5 = sig_new[ipos2:ipos2 + 6]
            p_t = (st0 + st1 + st2) / 3.0
            p_n = (sn0 + sn1 + sn2) / 3.0
            st0 -= p_t
            st1 -= p_t
            st2 -= p_t
            sn0 -= p_n
            sn1 -= p_n
            sn2 -= p_n
            sig_mises_test = np.sqrt(1.5 * (st0 ** 2 + st1 ** 2 + st2 ** 2) +
                                     3.0 * (st3 ** 2 + st4 ** 2 + st5 ** 2))
            sig_mises_new = np.sqrt(1.5 * (sn0 ** 2 + sn1 ** 2 + sn2 ** 2) +
                                    3.0 * (sn3 ** 2 + sn4 ** 2 + sn5 ** 2))
            if sig_mises_new > 0.0:
                T = p_n / sig_mises_new
            else:
                T = 1.0e12

            critical_strain = alpha * np.exp(-beta * T)

            if (sig_mises_test > sig_yield):
                peeq[ipos1] += (sig_mises_test - sig_yield) / G / 3.0

            csr[ipos1] = peeq[ipos1] / critical_strain


# update stresses and loads

@jit(nopython=True, cache=True, fastmath=True)
def update_stress_load(gp10, elNodes, nocoord, materialbyElement, sig_yield, du, sig, sig_update, sig_test_global, qin):
    u10 = np.empty(30, dtype=np.float64)  # displacements for the 10 tetrahedral nodes
    sig_test = np.empty(6, dtype=np.float64)
    bm0 = np.empty(10, dtype=np.float64)
    bm1 = np.empty(10, dtype=np.float64)
    bm2 = np.empty(10, dtype=np.float64)
    bm3 = np.empty(10, dtype=np.float64)
    bm4 = np.empty(10, dtype=np.float64)
    bm5 = np.empty(10, dtype=np.float64)
    bm6 = np.empty(10, dtype=np.float64)
    bm7 = np.empty(10, dtype=np.float64)
    bm8 = np.empty(10, dtype=np.float64)
    dmat = np.zeros((6, 6), dtype=np.float64)
    xlv0 = np.empty(10, dtype=np.float64)
    xlv1 = np.empty(10, dtype=np.float64)
    xlv2 = np.empty(10, dtype=np.float64)
    dshp0 = np.empty(10, dtype=np.float64)
    dshp1 = np.empty(10, dtype=np.float64)
    dshp2 = np.empty(10, dtype=np.float64)
    dshpg0 = np.empty(10, dtype=np.float64)
    dshpg1 = np.empty(10, dtype=np.float64)
    dshpg2 = np.empty(10, dtype=np.float64)

    hooke(0, materialbyElement, dmat)  # local material stiffness matrix

    for el, nodes in enumerate(elNodes):
        for i, nd in enumerate(nodes):
            co = nocoord[nd - 1]
            xlv0[i] = co[0]
            xlv1[i] = co[1]
            xlv2[i] = co[2]
            # xlv[i] = co
        elpos = 24 * el  # 4 integration points with 6 stress components each
        elv = np.zeros(30, dtype=np.float64)  # element load vector, 10 nodes with 3 load components each

        for index, nd in enumerate(nodes):
            n3 = 3 * (nd - 1)
            i3 = 3 * index
            u10[i3] = du[n3]
            u10[i3 + 1] = du[n3 + 1]
            u10[i3 + 2] = du[n3 + 2]
        for i in range(4):
            ip = gp10[i]
            ippos = elpos + 6 * i

            # ------------------------------------------------------------------------------------------------------------
            xi = ip[0]
            et = ip[1]
            ze = ip[2]

            # local derivatives of the shape functions: xi-derivative - source: Calculix, G Dhondt
            dshp0[0] = 1.0 - 4.0 * (1.0 - xi - et - ze)
            dshp0[1] = 4.0 * xi - 1.0
            dshp0[2] = 0.0
            dshp0[3] = 0.0
            dshp0[4] = 4.0 * (1.0 - 2.0 * xi - et - ze)
            dshp0[5] = 4.0 * et
            dshp0[6] = -4.0 * et
            dshp0[7] = -4.0 * ze
            dshp0[8] = 4.0 * ze
            dshp0[9] = 0.0

            # local derivatives of the shape functions: eta-derivative - source: Calculix, G Dhondt
            dshp1[0] = 1.0 - 4.0 * (1.0 - xi - et - ze)
            dshp1[1] = 0.0
            dshp1[2] = 4.0 * et - 1.0
            dshp1[3] = 0.0
            dshp1[4] = -4.0 * xi
            dshp1[5] = 4.0 * xi
            dshp1[6] = 4.0 * (1.0 - xi - 2.0 * et - ze)
            dshp1[7] = -4.0 * ze
            dshp1[8] = 0.0
            dshp1[9] = 4.0 * ze

            # local derivatives of the shape functions: zeta-derivative - source: Calculix, G Dhondt
            dshp2[0] = 1.0 - 4.0 * (1.0 - xi - et - ze)
            dshp2[1] = 0.0
            dshp2[2] = 0.0
            dshp2[3] = 4.0 * ze - 1.0
            dshp2[4] = -4.0 * xi
            dshp2[5] = 0.0
            dshp2[6] = -4.0 * et
            dshp2[7] = 4.0 * (1.0 - xi - et - 2.0 * ze)
            dshp2[8] = 4.0 * xi
            dshp2[9] = 4.0 * et

            # local derivative of the global coordinates
            xs00 = xs01 = xs02 = xs10 = xs11 = xs12 = xs20 = xs21 = xs22 = 0.0
            for kb in range(10):
                xlv0kb = xlv0[kb]
                xlv1kb = xlv1[kb]
                xlv2kb = xlv2[kb]
                dshp0kb = dshp0[kb]
                dshp1kb = dshp1[kb]
                dshp2kb = dshp2[kb]
                xs00 += xlv0kb * dshp0kb
                xs01 += xlv0kb * dshp1kb
                xs02 += xlv0kb * dshp2kb
                xs10 += xlv1kb * dshp0kb
                xs11 += xlv1kb * dshp1kb
                xs12 += xlv1kb * dshp2kb
                xs20 += xlv2kb * dshp0kb
                xs21 += xlv2kb * dshp1kb
                xs22 += xlv2kb * dshp2kb

            # Jacobian
            xsj = (xs00 * xs11 * xs22 -
                   xs00 * xs12 * xs21 +
                   xs02 * xs10 * xs21 -
                   xs02 * xs11 * xs20 +
                   xs01 * xs12 * xs20 -
                   xs01 * xs10 * xs22)

            # global derivative of the local coordinates
            xsi00 = (xs11 * xs22 - xs21 * xs12) / xsj
            xsi01 = (xs02 * xs21 - xs01 * xs22) / xsj
            xsi02 = (xs01 * xs12 - xs02 * xs11) / xsj
            xsi10 = (xs12 * xs20 - xs10 * xs22) / xsj
            xsi11 = (xs00 * xs22 - xs02 * xs20) / xsj
            xsi12 = (xs10 * xs02 - xs00 * xs12) / xsj
            xsi20 = (xs10 * xs21 - xs20 * xs11) / xsj
            xsi21 = (xs20 * xs01 - xs00 * xs21) / xsj
            xsi22 = (xs00 * xs11 - xs10 * xs01) / xsj

            # global derivatives of the shape functions
            for jb in range(10):
                dshpg0[jb] = xsi00 * dshp0[jb] + xsi10 * dshp1[jb] + xsi20 * dshp2[jb]
                dshpg1[jb] = xsi01 * dshp0[jb] + xsi11 * dshp1[jb] + xsi21 * dshp2[jb]
                dshpg2[jb] = xsi02 * dshp0[jb] + xsi12 * dshp1[jb] + xsi22 * dshp2[jb]

            # strain interpolation matrix bm1 .. bm8
            for ib in range(10):
                d00 = dshpg0[ib]
                d10 = dshpg1[ib]
                d20 = dshpg2[ib]
                bm0[ib] = d00
                bm1[ib] = d10
                bm2[ib] = d20
                bm3[ib] = d10
                bm4[ib] = d00
                bm5[ib] = d20
                bm6[ib] = d00
                bm7[ib] = d20
                bm8[ib] = d10

            # ------------------------------------------------------------------------------------------

            # strain
            eps = np.zeros(6, dtype=np.float64)
            for j in range(10):
                j3 = 3 * j
                eps[0] += bm0[j] * u10[j3]
                eps[1] += bm1[j] * u10[j3 + 1]
                eps[2] += bm2[j] * u10[j3 + 2]
                eps[3] += bm3[j] * u10[j3] + bm4[j] * u10[j3 + 1]
                eps[4] += bm5[j] * u10[j3] + bm6[j] * u10[j3 + 2]
                eps[5] += bm7[j] * u10[j3 + 1] + bm8[j] * u10[j3 + 2]

            # elastic test stress
            for j in range(6):
                tmp = sig[ippos + j]
                for k in range(6):
                    tmp += dmat[j, k] * eps[k]
                sig_test[j] = tmp

            sig_test_global[ippos:ippos + 6] = sig_test
            # return elastic stress to von Mises yield surface
            sxx, syy, szz, sxy, syz, szx = vmises_original_optimised(sig_test, sig_yield)
            sig_update[ippos:ippos + 6] = sxx, syy, szz, sxy, syz, szx

            # calculate element load vectors
            ipxsj = ip[3] * abs(xsj)
            for j in range(10):
                j3 = 3 * j
                elv[j3] += (bm0[j] * sxx + bm3[j] * sxy + bm5[j] * syz) * ipxsj
                elv[j3 + 1] += (bm1[j] * syy + bm4[j] * sxy + bm7[j] * szx) * ipxsj
                elv[j3 + 2] += (bm2[j] * szz + bm6[j] * syz + bm8[j] * szx) * ipxsj

        # add element load vectors to the global load vector
        for i in range(10):
            iglob = nodes[i] - 1
            iglob3 = 3 * iglob
            i3 = 3 * i
            for k in range(3):
                qin[iglob3 + k] += elv[i3 + k]

    return


@jit("types.UniTuple(float64, 6)(float64[::1], float64)", nopython=True, cache=True)
def vmises_original_optimised(sig_test, sig_yield):
    st0, st1, st2, st3, st4, st5 = sig_test
    p = (st0 + st1 + st2) / 3.0
    st0 -= p
    st1 -= p
    st2 -= p
    sig_mises = np.sqrt(1.5 * (st0 ** 2 + st1 ** 2 + st2 ** 2) +
                        3.0 * (st3 ** 2 + st4 ** 2 + st5 ** 2))

    if sig_yield > sig_mises:
        fac = 1.0
    else:
        fac = sig_yield / sig_mises

    st0 = fac * st0 + p
    st1 = fac * st1 + p
    st2 = fac * st2 + p
    st3 = fac * st3
    st4 = fac * st4
    st5 = fac * st5

    sig_update = st0, st1, st2, st3, st4, st5

    return sig_update


# map stresses to nodes
def mapStresses(elNodes, nocoord, sig, peeq, noce):
    # map maps corner node stresses to all tet10 nodes
    map = np.array([[1.0, 0.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0, 1.0],
                    [0.5, 0.5, 0.0, 0.0],
                    [0.0, 0.5, 0.5, 0.0],
                    [0.5, 0.0, 0.5, 0.0],
                    [0.5, 0.0, 0.0, 0.5],
                    [0.0, 0.5, 0.0, 0.5],
                    [0.0, 0.0, 0.5, 0.5]])

    expm = np.zeros((4, 4), dtype=np.float64)  # extrapolation matrix from Gauss points to corner nodes
    ipstress = np.zeros((4, 6), dtype=np.float64)  # Tet10 stresses by Gauss point (4 GP and 6 components)
    ippeeq = np.zeros((4, 1), dtype=np.float64)  # Tet10 peeq by Gauss point (4 GP and 1 component)

    ip10, ip6 = gaussPoints()

    tet10stress = np.zeros((len(nocoord), 6), dtype=np.float64)
    tet10peeq = np.zeros((len(nocoord)), dtype=np.float64)

    # map stresses in volumes to nodal points
    for el, nodes in enumerate(elNodes):
        elpos = 24 * el
        elposeq = 4 * el
        xl = np.array([nocoord[nd - 1] for nd in nodes]).T
        for index, ip in enumerate(ip10):
            xi = ip[0]
            et = ip[1]
            ze = ip[2]
            shp = shape4tet(xi, et, ze, xl)
            ippos = elpos + 6 * index
            ipposeq = elposeq + index
            ipstress[index] = sig[ippos:ippos + 6]  # ipstress (4x6): 6 stress components for 4 integration points
            ippeeq[index] = peeq[ipposeq]  # ippeeq (4x1): 1 strain component for 4 integration points
            for i in range(4):
                expm[index, i] = shp[i]
        expm_inv = np.linalg.inv(expm)
        npstress4 = np.dot(expm_inv, ipstress)  # npstress4 (4x6): for each corner node (4) all stress components (6)
        nppeeq4 = np.dot(expm_inv, ippeeq)  # nppeeq (4x1): for each corner node (4) one strain component (1)
        numnodes = np.array(
            [noce[nodes[n] - 1] for n in range(10)])  # numnodes = number of elements connected to node "nodes[n]-1"
        npstress10 = np.divide(np.dot(map, npstress4).T,
                               numnodes).T  # nodal point stress all nodes divided by number of connecting elements
        nppeeq10 = np.divide(np.dot(map, nppeeq4).T,
                             numnodes).T  # nodal point strain all nodes divided by number of connecting elements
        for index, nd in enumerate(nodes):
            tet10stress[nd - 1] += npstress10[index]

            tet10peeq[nd - 1] += nppeeq10[index]

    return tet10stress, tet10peeq


# fill resultobject with results
def pasteResults(doc, elNodes, nocoord, dis, tet10stress, tet10peeq):
    analysis = doc.getObject("Analysis")

    if analysis is None:
        prn_upd("No Analysis object. Please create one first")
        raise SystemExit()

    nn = len(nocoord)  # number of nodes

    resVol = analysis.addObject(ObjectsFem.makeResultMechanical(doc))[0]

    # VOLUME MESH START
    elements_tetra10 = {}
    mode_results_vol = {}

    nodes = range(1, nn + 1)
    coordinates = map(App.Vector, nocoord)
    displacements = map(App.Vector, np.array_split(dis, nn))
    volnodes = dict(zip(nodes, coordinates))
    mode_disp_vol = dict(zip(nodes, displacements))

    for index, elem in enumerate(elNodes):
        elements_tetra10[index + 1] = (
            elem[0], elem[2], elem[1], elem[3], elem[6], elem[5], elem[4], elem[7], elem[9], elem[8])

    mode_results_vol['disp'] = mode_disp_vol

    results = [mode_results_vol]

    mvol = {
        'Nodes': volnodes,
        'Seg2Elem': {},
        'Seg3Elem': {},
        'Tria3Elem': {},
        'Tria6Elem': {},
        'Quad4Elem': {},
        'Quad8Elem': {},
        'Tetra4Elem': {},
        'Tetra10Elem': elements_tetra10,
        'Hexa8Elem': {},
        'Hexa20Elem': {},
        'Penta6Elem': {},
        'Penta15Elem': {},
        'Results': results
    }

    meshvol = itf.make_femmesh(mvol)

    result_mesh_object_1 = ObjectsFem.makeMeshResult(doc, 'Result_Mesh_Volume')
    result_mesh_object_1.FemMesh = meshvol

    resVol.DisplacementVectors = [App.Vector(dis[3 * n], dis[3 * n + 1], dis[3 * n + 2]) for n in range(nn)]
    resVol.DisplacementLengths = [np.linalg.norm([dis[3 * n], dis[3 * n + 1], dis[3 * n + 2]]) for n in range(nn)]
    resVol.NodeStressXX = tet10stress.T[0].T.tolist()
    resVol.NodeStressYY = tet10stress.T[1].T.tolist()
    resVol.NodeStressZZ = tet10stress.T[2].T.tolist()
    resVol.NodeStressXY = tet10stress.T[3].T.tolist()
    resVol.NodeStressXZ = tet10stress.T[4].T.tolist()
    resVol.NodeStressYZ = tet10stress.T[5].T.tolist()
    resVol.Peeq = tet10peeq.T.tolist()

    resVol.Mesh = result_mesh_object_1
    resVol.NodeNumbers = [int(key) for key in resVol.Mesh.FemMesh.Nodes.keys()]

    resVol = itf.fill_femresult_mechanical(resVol, results)

    # VOLUME MESH FINISH

    # Add von Mises Stress and Plastic Strain Ratio to the results
    rt.add_von_mises(resVol)
    rt.add_principal_stress_std(resVol)

    doc.recompute()

    return resVol
