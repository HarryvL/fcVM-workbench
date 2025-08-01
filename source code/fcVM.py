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

import os
import sys
import time
import math
import dummyVM
import FemGui
import FreeCAD
import FreeCADGui
import ObjectsFem
import numpy as np
import Part as Part
import pyvista as pv
import FreeCAD as App
import FreeCADGui as Gui
from FreeCAD import Units
import scipy.sparse as scsp
from numba import jit, types
from numba.typed import Dict
from pyvista import CellType
import matplotlib.pyplot as plt
from femtools import membertools
from femmesh import meshsetsgetter
from femmesh import meshtools as mt
from scipy.sparse.linalg import eigsh, lobpcg
from scipy.linalg import eigh
from femresult import resulttools as rt
from feminout import importToolsFem as itf
from matplotlib.widgets import Button, TextBox
from matplotlib.ticker import FormatStrFormatter
from femtaskpanels import task_result_mechanical as trm

global mdir
mdir = os.path.dirname(dummyVM.file_path())
print("fcVM.py")

settings = {}
try:
    with open(os.path.join(mdir, 'fcVM.ini'), "r") as f:
        key = str(f.readline().strip()).split(" #")[0]
        settings[key] = int(f.readline().strip())
except FileNotFoundError:
    print("File fcVM.ini not found")

if settings["solver"] == 1:
    from sksparse.cholmod import cholesky
elif settings["solver"] == 2:
    from cholespy import CholeskySolverD, MatrixType
elif settings["solver"] == 3:
    from sksparse_minimal import SparseCholesky

global name
name = App.ActiveDocument.Label
file_path = os.path.join(mdir, "control files", name + '.inp')
macro_path = os.path.join(mdir, 'source code')
np.set_printoptions(precision=5, linewidth=300)


def prn_upd(*args):
    for obj in args:
        print(str(obj), end='')
    print('\n')
    Gui.updateGui()


def setUpAnalysis():
    doc = App.ActiveDocument

    return_code = 0

    print("hello from fcVM setUpAnalysis")

    mesh = None

    for obj in doc.Objects:
        if obj.Name[:7] == "FEMMesh":
            mesh = doc.getObject(obj.Name)

    if mesh is None:
        return_code = 1
    elif mesh.FemMesh.Nodes == {}:
        return_code = 2

    analysis = None

    for obj in doc.Objects:
        if obj.Name[:8] == "Analysis":
            analysis = doc.getObject(obj.Name)

    if analysis is None:
        return_code = 3

    # purge result objects
    if return_code == 0: rt.purge_results(analysis)

    doc.recompute()

    return doc, mesh, analysis, return_code


def setUpInput(doc, mesh, analysis):
    solver = doc.getObject("SolverCcxTools")
    member = membertools.AnalysisMember(analysis)

    if solver == None:
        FemGui.setActiveAnalysis(App.activeDocument().Analysis)
        if int(App.Version()[0])<1:
            FemGui.getActiveAnalysis().addObject(ObjectsFem.makeSolverCalculixCcxTools(App.ActiveDocument))
        else:
            FemGui.getActiveAnalysis().addObject(ObjectsFem.makeSolverCalculiXCcxTools(App.ActiveDocument))

        solver = doc.getObject("SolverCcxTools")

    # determine elements connected to a node using FC API
    fet = mt.get_femelement_table(mesh.FemMesh)
    # fet is dictionary: { elementid : [ nodeid, nodeid, ... , nodeid ] }
    net = mt.get_femnodes_ele_table(mesh.FemMesh.Nodes, fet)
    # net is dictionary: {nodeID : [[eleID, binary node position], [], ...], nodeID : [[], [], ...], ...}
    # node0 has binary node position 2^0 = 1, node1 = 2^1 = 2, ..., node10 = 2^10 = 1024

    # create nodal coordinate array nocoord for node number -> (x,y,z)
    ncv = list(mesh.FemMesh.Nodes.values())
    nocoord = np.asarray([[v.x, v.y, v.z] for v in ncv])  # nocoord[nodeIndex] = [x-coord, y-coord, z-coord]

    # get access to element sets: meshdatagetter.mat_geo_sets
    meshdatagetter = meshsetsgetter.MeshSetsGetter(
        analysis,
        solver,
        mesh,
        member)
    meshdatagetter.get_mesh_sets()

    if len(member.mats_linear) == 1:
        element_sets = [mesh.FemMesh.Volumes]
    else:
        element_sets = [es["FEMElements"] for es in member.mats_linear]

    # create connectivity array elNodes for mapping local node number -> global node number
    # elNodes = np.array(
    #     [mesh.FemMesh.getElementNodes(el) for el in mesh.FemMesh.Volumes])  # elNodes[elementIndex] = [node1,...,Node10]
    elNodes = np.array(
        [mesh.FemMesh.getElementNodes(el) for elset in element_sets for el in
         elset])  # elNodes[elementIndex] = [node1,...,Node10]

    matCon = {}  # BooleanFragment Primitive the material object refers to
    ppEl = {}  # BooleanFragment Primitive element El belongs to
    materialbyElement = []  # see further

    # prn_upd("Number of material objects: ", len(member.mats_linear))

    for indm, matobject in enumerate(member.mats_linear):
        E = float(App.Units.Quantity(matobject['Object'].Material['YoungsModulus']).getValueAs('MPa'))
        Nu = float(matobject['Object'].Material['PoissonRatio'])
        Density = float(App.Units.Quantity(matobject['Object'].Material['Density']).getValueAs('kg/mm^3'))
        # prn_upd("Material Object: ", matobject['Object'].Name, "   E= ", E, "   Nu= ", Nu, "   Density= ", Density)
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
                        bc = mesh.FemMesh.getNodesByVertex(ref)
                        for bcn in bc: bcnodes.append(bcn)
                    elif type(ref) == Part.Edge:
                        bc = mesh.FemMesh.getNodesByEdge(ref)
                        for bcn in bc: bcnodes.append(bcn)
                    elif type(ref) == Part.Face:
                        bc = mesh.FemMesh.getNodesByFace(
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
                try:
                    val = face[2][0].Value
                except:
                    val = face[2][0]
                fix[dof] = val
                fixdof[dof] = 0
        if not face[1][1]:
            for node in face[0]:
                dof = 3 * (node - 1) + 1
                try:
                    val = face[2][1].Value
                except:
                    val = face[2][1]
                fix[dof] = val
                fixdof[dof] = 0
        if not face[1][2]:
            for node in face[0]:
                dof = 3 * (node - 1) + 2
                try:
                    val = face[2][2].Value
                except:
                    val = face[2][2]
                fix[dof] = val
                fixdof[dof] = 0

    for dof in fix:
        if fix[dof] != 0.0:
            movdof[dof] = 1

    lf = [[0, 0, 0, 0, 0, 0]]  # load face nodes - signature for numba
    pr = [0.0]  # load face pressure - signature for numba
    lf_vertex = [[0]]
    pr_vertex = [[0.0, 0.0, 0.0]]
    lf_edge = [[0, 0, 0]]
    pr_edge = [[0.0, 0.0, 0.0]]
    lf_face = [[0, 0, 0, 0, 0, 0]]
    pr_face = [[0.0, 0.0, 0.0]]

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
                        for faceID in mesh.FemMesh.getFacesByFace(ref):  # face ID: ID of a 6-node face element
                            face_nodes = list(mesh.FemMesh.getElementNodes(faceID))  # 6-node element node numbers
                            lf.append(face_nodes)
                            if int(App.Version()[0]) < 1 and int(App.Version()[1]) < 22:
                                pr.append(sign * obj.Pressure)
                            else:
                                pr.append(sign * float(App.Units.Quantity(obj.Pressure.getValueAs('MPa'))))
                    else:
                        prn_upd("No Faces with Pressure Loads")

        if obj.isDerivedFrom('Fem::ConstraintForce'):
            if int(App.Version()[0]) < 1 and int(App.Version()[1]) < 22:
                F = obj.Force
            else:
                F = float(App.Units.Quantity(obj.Force.getValueAs('N')))
            d = obj.DirectionVector
            N = 0
            L = 0.0
            A = 0.0
            for part, boundaries in obj.References:
                for boundary in boundaries:
                    ref = part.Shape.getElement(boundary)
                    if type(ref) == Part.Vertex:
                        N += 1
                    elif type(ref) == Part.Edge:
                        L += ref.Length
                    else:
                        A += ref.Area

            for part, boundaries in obj.References:
                for boundary in boundaries:
                    ref = part.Shape.getElement(boundary)
                    if type(ref) == Part.Vertex:
                        dp = [F * d.x / N, F * d.y / N, F * d.z / N]
                        lf_vertex.append(list(mesh.FemMesh.getNodesByVertex(ref)))
                        pr_vertex.append(dp)
                    elif type(ref) == Part.Edge:
                        dl = [F * d.x / L, F * d.y / L, F * d.z / L]
                        for edgeID in mesh.FemMesh.getEdgesByEdge(ref):
                            lf_edge.append(list(mesh.FemMesh.getElementNodes(edgeID)))
                            pr_edge.append(dl)
                    elif type(ref) == Part.Face:
                        dp = [F * d.x / A, F * d.y / A, F * d.z / A]
                        for faceID in mesh.FemMesh.getFacesByFace(ref):
                            lf_face.append(list(mesh.FemMesh.getElementNodes(faceID)))
                            pr_face.append(dp)
                    else:
                        prn_upd("No Boundaries Found")

    loadfaces = np.array(lf)
    pressure = np.array(pr)
    loadvertices = np.array(lf_vertex)
    vertexloads = np.array(pr_vertex)
    loadedges = np.array(lf_edge)
    edgeloads = np.array(pr_edge)
    loadfaces_uni = np.array(lf_face)
    faceloads = np.array(pr_face)

    # re-order element nodes
    for el in elNodes:
        el[1], el[2] = el[2], el[1]
        el[4], el[6] = el[6], el[4]
        el[8], el[9] = el[9], el[8]

    return (elNodes, nocoord, fix, fixdof, movdof, materialbyElement, noce,
            loadfaces, pressure,
            loadvertices, vertexloads,
            loadedges, edgeloads,
            loadfaces_uni, faceloads)


# shape functions for a 4-node tetrahedron - only used for stress interpolation
@jit(nopython=True, cache=True)
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


@jit(nopython=True, cache=True, fastmath=True)
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

    # lmat = dshpg.T

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

    return xsj, dshpg


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


@jit(nopython=True, cache=True)
def shape2lin(xi, xle):
    shp = np.zeros((3), dtype=np.float64)
    dshp = np.zeros((3), dtype=np.float64)

    # shape functions
    shp[0] = - 0.5 * (1.0 - xi) * xi
    shp[1] = 0.5 * (1.0 + xi) * xi
    shp[2] = (1.0 + xi) * (1.0 - xi)

    # local derivatives of the shape functions: xi-derivative
    dshp[0] = xi - 0.5
    dshp[1] = xi + 0.5
    dshp[2] = - 2.0 * xi

    dx_dxi = xle[0][0] * dshp[0] + xle[0][1] * dshp[1] + xle[0][2] * dshp[2]
    dy_dxi = xle[1][0] * dshp[0] + xle[1][1] * dshp[1] + xle[1][2] * dshp[2]
    dz_dxi = xle[2][0] * dshp[0] + xle[2][1] * dshp[1] + xle[2][2] * dshp[2]

    xsj = np.sqrt(dx_dxi ** 2 + dy_dxi ** 2 + dz_dxi ** 2)

    return xsj, shp


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
    # Gaussian integration points and weights for 3-noded line
    gp2 = np.array([[-0.5773502691896257, 1.0], [0.5773502691896257, 1.0]])

    return gp10, gp6, gp2


# calculate the global stiffness matrix and load vector
# @jit(
#     "types.Tuple((float64[:],int64[:], int64[:], float64[:], float64[:]))(int64[:,:], float64[:,:], float64[:,:], DictType(int64,float64), int64[:,:], float64, float64, float64, float64[:])",
#     nopython=True, cache=True)
@jit(nopython=True, cache=True)
def calcGSM(elNodes, nocoord, materialbyElement, fix, grav_x, grav_y, grav_z, loadfaces, pressure,
            loadvertices, vertexloads, loadedges, edgeloads, loadfaces_uni, faceloads):
    gp10, gp6, gp2 = gaussPoints()
    ne = len(elNodes)  # number of volume elements
    # nn = len(nocoord[:, 0])  # number of degrees of freedom
    nn = len(nocoord)  # number of degrees of freedom

    # number of entries in the lower diagonal of the stiffness matrix: (dof**2 + dof) / 2
    ns = int((30 * 30 + 30) / 2 * ne)

    print(f"ne = {ne}")
    print(f"ns = {ns}\n")

    xle = np.zeros((3, 3), dtype=np.float64)  # coordinates of load line nodes
    xlf = np.zeros((3, 6), dtype=np.float64)  # coordinates of load face nodes
    xlv = np.zeros((10, 3), dtype=np.float64)  # coordinates of volume element nodes
    glv = np.zeros((3 * nn), dtype=np.float64)  # global load vector
    row = np.zeros(ns, dtype=np.int64)  # row indices of COO matrix
    col = np.zeros(ns, dtype=np.int64)  # column indices of COO matrix
    stm = np.zeros(ns, dtype=np.float64)  # stiffness values of COO matrix
    modf = np.zeros((3 * nn), dtype=np.float64)  # modification to the stiffness matrix for displacement BCs
    dof = np.zeros(30, dtype=np.int64)
    dmat = np.zeros((6, 6), dtype=np.float64)
    bmatV = np.zeros((6, 30), dtype=np.float64)
    x = np.zeros((4 * ne, 3), dtype=np.float64)

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

    for vertex in range(len(loadvertices) - 1):  # first vertex is a dummy signature for numba
        if len(loadvertices) == 1:
            break

        nd = loadvertices[vertex + 1][0]
        iglob = nd - 1
        iglob3 = 3 * iglob
        glv[iglob3:iglob3 + 3] += vertexloads[vertex + 1]

    for face in range(len(loadfaces_uni) - 1):  # first face is a dummy signature for numba
        if len(loadfaces_uni) == 1:
            break

        nda = loadfaces_uni[face + 1]  # node numbers of loaded face
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
            for i in range(len(loadfaces_uni[face] - 1)):
                nd = loadfaces_uni[face + 1][i]
                iglob = nd - 1
                iglob3 = 3 * iglob
                load = shp[nl] * faceloads[face + 1] * abs(xsj) * gp6[index][2]
                glv[iglob3:iglob3 + 3] += load
                nl += 1

    for edge in range(len(loadedges) - 1):  # first edge is a dummy signature for numba
        if len(loadedges) == 1:
            break

        nda = loadedges[edge + 1]  # node numbers of loaded edge
        for i in range(3):
            for j in range(3):
                nd = nda[j]
                xle[i][j] = nocoord[nd - 1][i]  # coordinates of loaded edge nodes
        # integrate element load vector
        for index in range(len(gp2)):
            xi = gp2[index][0]
            xsj, shp = shape2lin(xi, xle)
            nl = 0
            for i in range(len(loadedges[edge] - 1)):
                nd = loadedges[edge + 1][i]
                iglob = nd - 1
                iglob3 = 3 * iglob
                load = shp[nl] * edgeloads[edge + 1] * abs(xsj) * gp2[index][1]
                glv[iglob3:iglob3 + 3] += load
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
        for i, ip in enumerate(gp10):
            xi = ip[0]
            et = ip[1]
            ze = ip[2]
            shp = shp10tet(xi, et, ze)
            xsj, dshpg = dshp10tet(xi, et, ze, xlv, bmatV)
            esm += np.dot(bmatV.T, np.dot(dmat, bmatV)) * ip[3] * abs(xsj)
            gamma[0::3] += grav_x * density * shp * ip[3] * abs(xsj)
            gamma[1::3] += grav_y * density * shp * ip[3] * abs(xsj)
            gamma[2::3] += grav_z * density * shp * ip[3] * abs(xsj)
            V += xsj * ip[3]  # Element volume - not used
            x[4 * el + i] = np.dot(xlv.T, shp)  #

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

    return stm, row, col, glv, modf, V, loadsumx, loadsumy, loadsumz, ne, nn, x


@jit(nopython=True, cache=True)
def calcTSM(nstep, elNodes, nocoord, materialbyElement, fix, grav_x, grav_y, grav_z, loadfaces, pressure,
            loadvertices, vertexloads, loadedges, edgeloads, loadfaces_uni, faceloads, disp_new, du, sig_old, pgp,
            Et_E):
    gp10, gp6, gp2 = gaussPoints()
    ne = len(elNodes)  # number of volume elements
    nn = len(nocoord)  # number of degrees of freedom

    # normal analysis - number of entries in the lower diagonal of the stiffness matrix
    ns = int((30 * 30 + 30) / 2 * ne)

    # linear buckling analysis - number of entries in the full stiffness matrix
    nslb = int(30 * 30 * ne)

    # print(f"update tangent stiffness matrix")
    xle = np.zeros((3, 3), dtype=np.float64)  # coordinates of load line nodes
    xlf = np.zeros((3, 6), dtype=np.float64)  # coordinates of load face nodes
    xlv = np.zeros((10, 3), dtype=np.float64)  # coordinates of volume element nodes
    glv = np.zeros((3 * nn), dtype=np.float64)  # global load vector
    if (float(nstep) > 1.0):
        row = np.zeros(ns, dtype=np.int64)  # row indices of COO matrix
        col = np.zeros(ns, dtype=np.int64)  # column indices of COO matrix
    else:
        row = np.zeros(nslb, dtype=np.int64)  # row indices of COO matrix
        col = np.zeros(nslb, dtype=np.int64)  # column indices of COO matrix
    stm = np.zeros(ns, dtype=np.float64)  # stiffness values of COO matrix
    stms = np.zeros(nslb, dtype=np.float64)  # stiffness values of COO matrix
    stmg = np.zeros(nslb, dtype=np.float64)  # stiffness values of COO matrix
    modf = np.zeros((3 * nn), dtype=np.float64)  # modification to the stiffness matrix for displacement BCs
    modfs = np.zeros((3 * nn), dtype=np.float64)  # modification to the stiffness matrix for displacement BCs
    modfg = np.zeros((3 * nn), dtype=np.float64)  # modification to the stiffness matrix for displacement BCs
    dof = np.zeros(30, dtype=np.int64)
    dmat = np.zeros((6, 6), dtype=np.float64)
    pmat = np.zeros((6, 6), dtype=np.float64)
    bmatV = np.zeros((6, 30), dtype=np.float64)
    x = np.zeros((4 * ne, 3), dtype=np.float64)

    #   calculate element load vectors for pressure and add to global vector

    for face in range(len(pressure) - 1):  # first pressure value is a dummy signature for numba
        if len(pressure) == 1:
            break

        nda = loadfaces[face + 1]  # node numbers of loaded face
        for i in range(3):
            for j in range(6):
                nd = nda[j]
                idof = 3 * (nd - 1) + i
                xlf[i][j] = nocoord[nd - 1][i] + disp_new[
                    idof]  # coordinates of loaded face nodes - pressure on stretched surface

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

    for vertex in range(len(loadvertices) - 1):  # first vertex is a dummy signature for numba
        if len(loadvertices) == 1:
            break

        nd = loadvertices[vertex + 1][0]
        iglob = nd - 1
        iglob3 = 3 * iglob
        glv[iglob3:iglob3 + 3] += vertexloads[vertex + 1]

    for face in range(len(loadfaces_uni) - 1):  # first face is a dummy signature for numba
        if len(loadfaces_uni) == 1:
            break

        nda = loadfaces_uni[face + 1]  # node numbers of loaded face
        for i in range(3):
            for j in range(6):
                nd = nda[j]
                xlf[i][j] = nocoord[nd - 1][i]  # coordinates of loaded face nodes - load on original surface

        # integrate element load vector
        for index in range(len(gp6)):
            xi = gp6[index][0]
            et = gp6[index][1]
            xsj, shp, bmatS, xx, xt, xp = shape6tri(xi, et, xlf)
            nl = 0
            for i in range(len(loadfaces_uni[face] - 1)):
                nd = loadfaces_uni[face + 1][i]
                iglob = nd - 1
                iglob3 = 3 * iglob
                load = shp[nl] * faceloads[face + 1] * abs(xsj) * gp6[index][2]
                glv[iglob3:iglob3 + 3] += load
                nl += 1

    for edge in range(len(loadedges) - 1):  # first edge is a dummy signature for numba
        if len(loadedges) == 1:
            break

        nda = loadedges[edge + 1]  # node numbers of loaded edge
        for i in range(3):
            for j in range(3):
                nd = nda[j]
                xle[i][j] = nocoord[nd - 1][i]  # coordinates of loaded edge nodes - load on original edge
        # integrate element load vector
        for index in range(len(gp2)):
            xi = gp2[index][0]
            xsj, shp = shape2lin(xi, xle)
            nl = 0
            for i in range(len(loadedges[edge] - 1)):
                nd = loadedges[edge + 1][i]
                iglob = nd - 1
                iglob3 = 3 * iglob
                load = shp[nl] * edgeloads[edge + 1] * abs(xsj) * gp2[index][1]
                glv[iglob3:iglob3 + 3] += load
                nl += 1

    # for each volume element calculate the element stiffness matrix
    # and gravity load vector and add to global matrix and vector

    pos = 0

    V = 0.0

    hooke(0, materialbyElement, dmat)
    density = materialbyElement[0][2]
    E = float(materialbyElement[0][0])
    nu = float(materialbyElement[0][1])
    G = E / (1.0 + nu) / 2.0
    if Et_E > 0.95: Et_E = 0.95
    Et = Et_E * E
    H = Et / (1.0 - Et_E)

    for el, nodes in enumerate(elNodes):
        esm = np.zeros((30, 30), dtype=np.float64)
        nsm = np.zeros((30, 30), dtype=np.float64)
        gamma = np.zeros((30), dtype=np.float64)
        elpos = 24 * el
        # set up nodal values for this element
        for i in range(3):
            for j in range(10):
                nd = nodes[j]
                idof = 3 * (nd - 1) + i
                xlv[j][i] = nocoord[nd - 1][i] + disp_new[
                    idof]  # + du[idof] # stretched coordinates (WARNING: gravity loads change)
                # TODO: gamma should be based on original coordinates and esm on updated coordinates

        # integrate element matrix
        for i, ip in enumerate(gp10):
            ip4 = 4 * el + i
            ip24 = 24 * el + 6 * i
            sig0 = sig_old[ip24:ip24 + 6]
            sig = np.array([[sig0[0], sig0[3], sig0[4]],
                            [sig0[3], sig0[1], sig0[5]],
                            [sig0[4], sig0[5], sig0[2]]])
            xi = ip[0]
            et = ip[1]
            ze = ip[2]
            shp = shp10tet(xi, et, ze)
            xsj, dshpg = dshp10tet(xi, et, ze, xlv, bmatV)
            if pgp[ip4]:
                st0, st1, st2, st3, st4, st5 = sig0
                p = (st0 + st1 + st2) / 3.0
                st0 -= p
                st1 -= p
                st2 -= p
                svm = np.sqrt(1.5 * (st0 ** 2 + st1 ** 2 + st2 ** 2) +
                              3.0 * (st3 ** 2 + st4 ** 2 + st5 ** 2))
                s = np.array([st0, st1, st2, st3, st4, st5])
                if svm == 0.0: svm = 1.0  # prevent division error in pmat

                fac = 3.0 * G / (1.0 + H / 3.0 / G) / svm ** 2
                for i1 in range(6):
                    for i2 in range(6):
                        pmat[i1][i2] = fac * s[i1] * s[i2]
                esm += np.dot(bmatV.T, np.dot(dmat - pmat, bmatV)) * ip[3] * abs(xsj)
            else:
                esm += np.dot(bmatV.T, np.dot(dmat, bmatV)) * ip[3] * abs(xsj)

            if float(nstep) == 1.0:
                GM = np.kron(dshpg, np.eye(3))
                SM = np.kron(sig, np.eye(3))
                dnsm = np.dot(GM.T, np.dot(SM, GM)) * ip[3] * abs(xsj)
                nsm += dnsm
            gamma[0::3] += grav_x * density * shp * ip[3] * abs(xsj)
            gamma[1::3] += grav_y * density * shp * ip[3] * abs(xsj)
            gamma[2::3] += grav_z * density * shp * ip[3] * abs(xsj)
            V += xsj * ip[3]  # Element volume - not used
            x[4 * el + i] = np.dot(xlv.T, shp)

        for i in range(10):
            nd = nodes[i] - 1
            glv[3 * nd] += gamma[3 * i]
            glv[3 * nd + 1] += gamma[3 * i + 1]
            glv[3 * nd + 2] += gamma[3 * i + 2]
            for j in range(3):
                dof[3 * i + j] = 3 * nd + j

        if (float(nstep) > 1.0):
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
                            modf[dofj] -= (esm[i][j] + nsm[i][j]) * fix[dofi]
                else:
                    for j in range(i + 1):
                        dofj = dof[j]
                        if dofj in fix:
                            modf[dofi] -= (esm[i][j] + nsm[i][j]) * fix[dofj]
                        else:
                            if dofi > dofj:
                                row[pos] = dofi
                                col[pos] = dofj
                            else:
                                row[pos] = dofj
                                col[pos] = dofi
                            stm[pos] = esm[i][j] + nsm[i][j]
                            # stm[pos] = esm[i][j]
                            stms[pos] = esm[i][j]
                            stmg[pos] = nsm[i][j]
                            pos += 1
        else:
            for i in range(30):
                dofi = dof[i]
                for j in range(30):
                    dofj = dof[j]
                    row[pos] = dofi
                    col[pos] = dofj
                    stms[pos] = esm[i][j]
                    stmg[pos] = nsm[i][j]
                    if i == j and dofi in fix:
                        stms[pos] *= 100.0
                    pos += 1

    row = row[:pos]
    col = col[:pos]
    stm = stm[:pos]
    stms = stms[:pos]
    stmg = stmg[:pos]

    loadsumx = 0.0
    loadsumy = 0.0
    loadsumz = 0.0
    for node in range(nn):
        dof = 3 * node
        loadsumx += glv[dof]
        loadsumy += glv[dof + 1]
        loadsumz += glv[dof + 2]

    return stm, stms, stmg, row, col, glv, modf


# calculate load-deflection curve
def calcDisp(elNodes, nocoord, fixdof, movdof, modf, materialbyElement, stm, row, col,
             glv, nstep, iterat_max, error_max, relax, scale_re, scale_up, scale_dn, sig_yield_inp, disp_output,
             ultimate_strain, fcVM_window, Et_E, target_LF, x, noce, fix, grav_x, grav_y, grav_z, loadfaces, pload,
             loadvertices, vertexloads, loadedges, edgeloads, loadfaces_uni, faceloads, gnl, maxImp, ev1, ev2):
    if gnl == "GNLY":
        LD = True
        relax = 1.0
        disp_output = "total"
        ev1 = float(ev1)
        ev2 = float(ev2)
        maxImp = float(maxImp)
        scale_up = 1.1

    else:
        LD = False

    ndof = len(glv)  # number of degrees of freedom
    nelem = len(elNodes)  # number of elements

    nocoord_old = np.copy(nocoord)

    pr_u = fcVM_window.progressBar.setValue
    st_u = fcVM_window.Step.setText
    LF_u = fcVM_window.Load_Factor.setText
    PQ_u = fcVM_window.PEEQ.setText
    CR_u = fcVM_window.CSR.setText

    t0 = time.perf_counter()
    gsm = scsp.csc_matrix((stm, (row, col)), shape=(ndof, ndof))  # construct sparse global stiffness matrix
    t1 = time.perf_counter()
    # prn_upd("construct sparse global stiffness matrix: {:.2e} s".format((t1 - t0)))

    qnorm = np.linalg.norm(glv)
    if qnorm < 1.0: qnorm = 1.0

    # Cholesky decomposition of the global stiffness matrix and elastic solution using Cholmod
    t0 = time.perf_counter()
    if settings["solver"] == 1:
        factor = cholesky(gsm)
    elif settings["solver"] == 2:
        solver = CholeskySolverD(gsm.shape[0], gsm.indptr, gsm.indices, gsm.data, MatrixType.CSC)
    elif settings["solver"] == 3:
        sparse_cholesky = SparseCholesky(gsm)

    t1 = time.perf_counter()
    f = fixdof * glv + modf
    if settings["solver"] == 1:
        ue = factor(f)  # elastic solution
    elif settings["solver"] == 2:
        ue = np.empty(gsm.shape[0])
        solver.solve(f, ue)
    elif settings["solver"] == 3:
        ue = sparse_cholesky.solve_A(f)

    disp_el = ue.copy()

    t2 = time.perf_counter()
    # prn_upd("sparse Cholesky decomposition: {:.2e} s, elastic solution: {:.2e} s".format((t1 - t0), (t2 - t1)))

    # initiate analysis
    dl0 = 1.0 / nstep
    dl = dl0
    du = dl * ue

    sig_new = np.zeros(24 * nelem, dtype=np.float64)  # stress in Tet10
    sig_old = np.zeros(24 * nelem, dtype=np.float64)
    sig_yield = np.full(4 * nelem, sig_yield_inp, dtype=np.float64)  # yield stress in Tet10
    sig_test = np.zeros(24 * nelem, dtype=np.float64)
    peeq = np.zeros(4 * nelem, dtype=np.float64)  # equivalent plastic strain in Tet10
    triax = np.zeros(4 * nelem, dtype=np.float64)  # triaxiality in Tet10
    pressure = np.zeros(4 * nelem, dtype=np.float64)  # pressure in Tet10
    sigmises = np.zeros(4 * nelem, dtype=np.float64)  # von Mises stress in Tet10
    pgp = np.full(4 * nelem, False, dtype=bool)
    ecr = np.zeros(4 * nelem, dtype=np.float64)  # critical plastic strain in Tet10
    csr = np.zeros(4 * nelem, dtype=np.float64)  # critical strain ratio in Tet10
    disp_new = np.zeros(ndof, dtype=np.float64)  # elastic displacement results
    disp_old = np.zeros(ndof, dtype=np.float64)  # displacement results
    lbd = np.zeros(1, dtype=np.float64)  # load level
    rfl = np.zeros(1, dtype=np.float64)  # reaction force level (for displacement control)
    eigenval = np.zeros(1, dtype=np.float64)
    eigenvec = np.zeros((1, 1), dtype=np.float64)
    imper = np.zeros(ndof, dtype=np.float64)

    gp10, gp6, gp2 = gaussPoints()

    # determine elastic reaction force on moving boundary
    if max(movdof) == 1:
        qelastic = np.zeros(3 * len(nocoord), dtype=np.float64)
        update_stress_load(gp10, elNodes, nocoord, materialbyElement, sig_yield, disp_new, ue, sig_old, sig_new,
                           sig_test, qelastic, Et_E, LD, pgp)

        qelastic *= movdof
        qelnorm = np.linalg.norm(qelastic)
        qnorm = qelnorm
        sig_new = np.zeros(24 * nelem, dtype=np.float64)  # reset sig_new to zero

    step = -1
    cnt = True
    fail = False
    dof_max = np.argmax(np.abs(ue))

    un = [0.]
    csrplot = [0.]
    crip = [0]
    pplot = [0.]
    svmplot = [0.]
    triaxplot = [0.]
    peeqplot = [0.]
    peeqmax = [0.]
    ecrplot = [0.]
    lout = [0.]

    update_stress_load(gp10, elNodes, nocoord, materialbyElement, 1.0e6 * sig_yield, np.zeros(ndof, dtype=np.float64),
                       ue, sig_old, sig_new, sig_test, np.zeros(3 * len(nocoord), dtype=np.float64), Et_E, False,
                       pgp)

    # perform an elastic buckling analysis
    if LD and not (float(nstep) > 1.0 and maxImp == 0.0):
        # linear buckling analysis
        stm, stms, stmg, row, col, glv, modf = calcTSM(1.0, elNodes, nocoord, materialbyElement, fix, grav_x,
                                                       grav_y, grav_z, loadfaces, pload, loadvertices, vertexloads,
                                                       loadedges, edgeloads,
                                                       loadfaces_uni, faceloads, np.zeros(ndof, dtype=np.float64),
                                                       np.zeros(ndof, dtype=np.float64),
                                                       sig_new, pgp, Et_E)

        K = scsp.csc_matrix((stms, (row, col)), shape=(ndof, ndof))
        G = - scsp.csc_matrix((stmg, (row, col)), shape=(ndof, ndof))

        eigenval, eigenvec = eigsh(K, k=2, M=G, sigma=0.1, which='LM', mode='buckling')

        print("buckling load factors: ", eigenval)

    if float(nstep) == 1.0:  # elastic and linear buckling analysis
        step = 0
        out_disp = 1
        disp_new = ue
        lbd = np.append(lbd, 1.0)
        rfl = np.append(rfl, 1.0)
        un.append(np.max(np.abs(disp_new)))
        cnt = False
    elif LD and maxImp != 0.0:  # apply imperfection and reset analysis

        ua = ev1 / (ev1 + ev2) * eigenvec[:, 0] + ev2 / (ev1 + ev2) * eigenvec[:, 1]
        ub = ev1 / (ev1 + ev2) * eigenvec[:, 0] - ev2 / (ev1 + ev2) * eigenvec[:, 1]
        ma = np.max(np.abs(ua))
        mb = np.max(np.abs(ub))

        if ma > mb:
            imax = np.argmax(np.abs(ua))
            imper = maxImp / ma * np.sign(ua[imax]) * ua
        else:
            imax = np.argmax(np.abs(ub))
            imper = maxImp / mb * np.sign(ub[imax]) * ub

        # print("maximum imperfection: ", np.max(imper))

        nocoord += imper.reshape(-1, 3)

        stm, row, col, glv, modf, _, _, _, _, _, _, _ = calcGSM(elNodes, nocoord,
                                                                materialbyElement,
                                                                fix, grav_x, grav_y,
                                                                grav_z,
                                                                loadfaces, pload,
                                                                loadvertices,
                                                                vertexloads,
                                                                loadedges,
                                                                edgeloads,
                                                                loadfaces_uni,
                                                                faceloads)

        t0 = time.perf_counter()
        gsm = scsp.csc_matrix((stm, (row, col)), shape=(ndof, ndof))  # construct sparse global stiffness matrix
        t1 = time.perf_counter()
        # prn_upd("update sparse global stiffness matrix for imperfections: {:.2e} s".format((t1 - t0)))
        qnorm = np.linalg.norm(glv)
        if qnorm < 1.0: qnorm = 1.0

        # Cholesky decomposition of the global stiffness matrix and elastic solution using Cholmod
        t0 = time.perf_counter()
        if settings["solver"] == 1:
            factor = cholesky(gsm)
        elif settings["solver"] == 2:
            solver = CholeskySolverD(gsm.shape[0], gsm.indptr, gsm.indices, gsm.data, MatrixType.CSC)
        elif settings["solver"] == 3:
            sparse_cholesky = SparseCholesky(gsm)

        t1 = time.perf_counter()
        f = fixdof * glv + modf
        if settings["solver"] == 1:
            ue = factor(f)  # elastic solution
        elif settings["solver"] == 2:
            ue = np.empty(gsm.shape[0])
            solver.solve(f, ue)
        elif settings["solver"] == 3:
            ue = sparse_cholesky.solve_A(f)

        disp_el = ue.copy()

        t2 = time.perf_counter()
        # prn_upd("sparse Cholesky decomposition of update stiffness matrix: {:.2e} s, elastic solution: {:.2e} s".format(
        #     (t1 - t0), (t2 - t1)))

        # initiate analysis
        dl0 = 1.0 / nstep  # nstep == 1 execute an elastic analysis
        dl = dl0
        du = dl * ue

        sig_old = np.zeros(24 * nelem, dtype=np.float64)
        sig_test = np.zeros(24 * nelem, dtype=np.float64)
        disp_new = np.zeros(ndof, dtype=np.float64)  # elastic displacement results
        disp_old = np.zeros(ndof, dtype=np.float64)  # displacement results
        lbd = np.zeros(1, dtype=np.float64)  # load level

    factor_time_tot = 0.0
    factor_dec_tot = 0.0
    iterat_tot = 0
    mrr = False
    sig_new = np.zeros(24 * nelem, dtype=np.float64)  # stress in Tet10
    pgp = np.full(4 * nelem, False, dtype=bool)

    while cnt:
        cnt = False
        iRiks = True
        pstep = 0
        pr_u(0)
        while pstep < nstep and not mrr:
            step += 1
            pstep += 1
            st_u(str(pstep))
            restart = 0

            prn_upd("Step: {}".format(step))
            a = du.copy()  # a: Riks control vector
            if iRiks:
                sig_old = sig_new.copy()
                lbd = np.append(lbd, lbd[step] + dl)  # lbd: load level
            else:
                lbd[step + 1] = lbd[step] + dl

            # update stresses and loads
            qin = np.zeros(3 * len(nocoord), dtype=np.float64)
            update_stress_load(gp10, elNodes, nocoord, materialbyElement, sig_yield, disp_new, du, sig_old,
                               sig_new, sig_test, qin, Et_E, LD, pgp)

            # calculate residual load vector
            fex = fixdof * lbd[step + 1] * glv
            fin = fixdof * qin
            r = fex - fin
            # print("at start of step")
            # print("fex: ", np.linalg.norm(fex))
            # print("fin: ", np.linalg.norm(fin))
            # print("r: ", np.linalg.norm(r))
            # print("sig_max: ", np.max(sig_new))

            rnorm = np.linalg.norm(r)

            # out-of-balance error
            error = rnorm / qnorm

            iterat = 0
            prn_upd("Iteration: {}, Error: {:.2e}".format(iterat, error))

            while error > error_max and not mrr:

                iterat += 1
                iterat_tot += 1

                if LD and (iterat == 1 or np.any(pgp)):

                    stm, stms, stmg, row, col, glv, modf = calcTSM(nstep, elNodes, nocoord, materialbyElement, fix,
                                                                   grav_x, grav_y, grav_z,
                                                                   loadfaces, pload, loadvertices, vertexloads,
                                                                   loadedges,
                                                                   edgeloads,
                                                                   loadfaces_uni, faceloads, disp_new, du, sig_old,
                                                                   pgp,
                                                                   Et_E)

                    t0 = time.perf_counter()
                    tsm = scsp.csc_matrix((stm, (row, col)),
                                          shape=(ndof, ndof))  # construct sparse global stiffness matrix
                    t1 = time.perf_counter()
                    # prn_upd("construct sparse tangent stiffness matrix: {:.2e} s".format((t1 - t0)))
                    try:
                        if settings["solver"] == 1:
                            factor = cholesky(tsm)
                        elif settings["solver"] == 2:
                            solver = CholeskySolverD(tsm.shape[0], tsm.indptr, tsm.indices, tsm.data, MatrixType.CSC)
                        elif settings["solver"] == 3:
                            sparse_cholesky = SparseCholesky(tsm)

                        t2 = time.perf_counter()
                        # prn_upd("decompose tangent stiffness matrix: {:.2e} s".format((t2 - t1)))
                        factor_dec_tot += t2 - t1

                    except:
                        print("singular stiffness matrix")
                        pass
                    fe = fixdof * glv + modf
                    if settings["solver"] == 1:
                        ue = factor(fe)  # elastic solution
                    elif settings["solver"] == 2:
                        ue = np.empty(tsm.shape[0])
                        solver.solve(fe, ue)
                    elif settings["solver"] == 3:
                        ue = sparse_cholesky.solve_A(fe)

                    t3 = time.perf_counter()

                    # prn_upd("solve for ue: {:.2e} s".format((t3 - t2)))

                    a = ue.copy()
                    a *= np.linalg.norm(du) / np.linalg.norm(a)

                f = relax * r
                t0 = time.perf_counter()
                if settings["solver"] == 1:
                    due = factor(f)
                elif settings["solver"] == 2:
                    due = np.empty(gsm.shape[0])
                    solver.solve(f, due)
                elif settings["solver"] == 3:
                    due = sparse_cholesky.solve_A(f)

                t1 = time.perf_counter()
                # prn_upd("solve for due: {:.2e} s".format((t1 - t0)))

                factor_time_tot += t1 - t0

                # Riks control correction to load level increment
                if iRiks:
                    dl = -np.dot(a, due) / np.dot(a, ue)
                    lbd[step + 1] += dl
                    aa = np.linalg.norm(a)
                    # print("||a||: ", aa)
                else:
                    dl = 0.0
                    # du += due + dl * ue

                # Riks control correction to displacement increment
                du += due + dl * ue

                uu = np.linalg.norm(du)
                # print("||du1||: ", uu)

                # scale back increment

                sf = min(aa / uu, 1.0)

                lbd[step + 1] = lbd[step] + sf * (lbd[step + 1] - lbd[step])
                du *= sf

                # uu = np.linalg.norm(du)
                # print("||du2||: ", uu)

                # update stresses and loads

                qin = np.zeros(3 * len(nocoord), dtype=np.float64)
                update_stress_load(gp10, elNodes, nocoord, materialbyElement, sig_yield, disp_new, du, sig_old,
                                   sig_new, sig_test, qin, Et_E, LD, pgp)

                # calculate out of balance error
                r = fixdof * (lbd[step + 1] * glv - qin)
                rnorm = np.linalg.norm(r)
                error = rnorm / qnorm
                # print("at end of iteration")
                # print("fex: ", np.linalg.norm(fixdof * lbd[step + 1] * glv))
                # print("fin: ", np.linalg.norm(fixdof * qin))
                # print("r: ", np.linalg.norm(r))
                # print("sig_max: ", np.max(sig_new))

                prn_upd("Iteration: {}, Error: {:.2e}".format(iterat, error))

                if iterat > iterat_max:
                    print(f"RESTART # {restart + 1}")
                    # scale down
                    if restart > 3:
                        print("MAXIMUM RESTARTS REACHED")
                        fail = False
                        step -= 1
                        lbd = lbd[:-1]  # remove the last element of lbd
                        mrr = True
                    restart += 1
                    if step > 0 and not mrr:
                        dl = (lbd[step] - lbd[step - 1]) / scale_re / restart
                        du = (disp_new - disp_old) / scale_re / restart
                    elif not mrr:
                        # for first step only
                        dl = dl0 / scale_re / restart
                        du = dl * ue / scale_re / restart
                    lbd[step + 1] = lbd[step] + dl

                    if not mrr:
                        qin = np.zeros(3 * len(nocoord), dtype=np.float64)
                        update_stress_load(gp10, elNodes, nocoord, materialbyElement, sig_yield, disp_new, du,
                                           sig_old, sig_new, sig_test, qin, Et_E, LD, pgp)
                        r = fixdof * (lbd[step + 1] * (glv + modf) - qin)
                        rnorm = np.linalg.norm(r)
                        error = rnorm / qnorm

                        iterat = 0

            if abs(target_LF - lbd[step]) < abs(lbd[step + 1] - lbd[step]) and iRiks:
                print("REACHED TARGET LOAD")
                fac = (target_LF - lbd[step]) / (lbd[step + 1] - lbd[step])
                du = fac * du
                sig_new = sig_old + fac * (sig_new - sig_old)
                sig_test = sig_old + fac * (sig_test - sig_old)
                lbd[step + 1] = target_LF
                disp_new += du
                disp_new_node = np.array(
                    [disp_new[3 * i] ** 2 + disp_new[3 * i + 1] ** 2 + disp_new[3 * i + 2] ** 2 for i in
                     range((ndof - 1) // 3)])
                un.append(np.sqrt(max(disp_new_node)))
                update_PEEQ_CSR(nelem, materialbyElement, sig_test, sig_new, sig_yield, ultimate_strain, peeq, csr,
                                triax,
                                pressure, sigmises, ecr, Et_E)
                maxloc = np.argmax(csr)
                csrplot.append(np.max(csr))
                crip.append(maxloc)
                pplot.append(pressure[maxloc])
                svmplot.append(sigmises[maxloc])
                triaxplot.append(triax[maxloc])
                ecrplot.append(ecr[maxloc])
                peeqplot.append(peeq[maxloc])
                peeqmax.append(np.max(peeq))
                break

                step -= 1
                pstep -= 1

            elif not mrr:
                # update results at end of converged load step
                pr_u(int(100 * (pstep + 1) / nstep))
                LF_u(str(round(lbd[step + 1], 3)))
                disp_old = disp_new.copy()
                disp_new += du
                dl = lbd[step + 1] - lbd[step]
                if max(movdof) == 1:
                    rfl = np.append(rfl, np.sum(movdof * qin))

                # TEST
                # dl = 0.0
                # du = np.zeros(ndof, dtype=np.float64)
                # TEST

                if iterat > 10:
                    # scale down
                    dl /= scale_dn
                    du /= scale_dn
                if iterat < 5:
                    # scale up
                    dl *= scale_up
                    du *= scale_up

                disp_new_node = np.array(
                    [disp_new[3 * i] ** 2 + disp_new[3 * i + 1] ** 2 + disp_new[3 * i + 2] ** 2 for i in
                     range((ndof - 1) // 3)])
                un.append(np.sqrt(max(disp_new_node)))
                update_PEEQ_CSR(nelem, materialbyElement, sig_test, sig_new, sig_yield, ultimate_strain, peeq, csr,
                                triax,
                                pressure, sigmises, ecr, Et_E)
                maxloc = np.argmax(csr)
                csrplot.append(np.max(csr))
                crip.append(maxloc)
                pplot.append(pressure[maxloc])
                svmplot.append(sigmises[maxloc])
                triaxplot.append(triax[maxloc])
                ecrplot.append(ecr[maxloc])
                peeqplot.append(peeq[maxloc])
                peeqmax.append(np.max(peeq))

                PQ_u(str(round(max(peeq), 3)))
                CR_u(str(round(max(csr), 3)))

                if not iRiks: break

        if max(movdof) == 1:
            lout = rfl
        else:
            lout = lbd

        print(
            '{0: >8}{1: >10}{2: >10}{3: >10}{4: >10}{5: >10}{6: >10}{7: >10}{8: >10}{9: >10}{10: >10}{11: >10}'.format(
                "ip_max",
                "x-coor",
                "y-coor",
                "z-coor",
                "load",
                "un",
                "peeq",
                "pressure",
                "svmises",
                "triax",
                "eps_cr",
                "csr_max"))
        for i in range(len(crip)):
            print(
                '{0: 8d}{1: >10.2e}{2: >10.2e}{3: >10.2e}{4: >10.2e}{5: >10.2e}{6: >10.2e}{7: >10.2e}{8: >10.2e}{9: >10.2e}{10: >10.2e}{11: >10.2e}'.format(
                    crip[i],
                    x[crip[i]][0],
                    x[crip[i]][1],
                    x[crip[i]][2],
                    lout[i],
                    un[i],
                    peeqplot[i],
                    pplot[i],
                    svmplot[i],
                    triaxplot[i],
                    ecrplot[i],
                    csrplot[i]))
        csr_non_zero = np.nonzero(csrplot)
        if len(csr_non_zero[0]) != 0:
            el_limit = csr_non_zero[0][0] - 1
        else:
            el_limit = 0
        csr_limit = np.argwhere(np.asarray(csrplot) > 1.0)
        peeq_limit = np.argwhere(np.asarray(peeqmax) > ultimate_strain)

        if fcVM_window.csrRbtn.isChecked():
            if len(csr_limit) != 0:
                ul_limit = csr_limit[0][0] - 1
            else:
                ul_limit = 0
        else:
            if len(peeq_limit) != 0:
                ul_limit = peeq_limit[0][0] - 1
            else:
                ul_limit = 0

        averaged = False
        cnt, dl, du, target_LF = plot(fcVM_window, averaged, el_limit, ul_limit, un, lout, csrplot, peeqmax, dl, du,
                                      target_LF,
                                      nstep, ue, ultimate_strain, disp_new, disp_old, elNodes, nocoord, sig_new, peeq,
                                      sigmises, csr, noce, sig_yield_inp)

    prn_upd("total time decomposing tangent stiffness matrix: {}".format(factor_dec_tot))
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
        return disp_new, disp_el, eigenval, eigenvec, sig_new, peeq, sigmises, csr, lout, un, crip, peeqplot, pplot, svmplot, triaxplot, ecrplot, csrplot, fail, nocoord_old
    else:
        return disp_new - disp_old, disp_el, eigenval, eigenvec, sig_new, peeq, sigmises, csr, lout, un, crip, peeqplot, pplot, svmplot, triaxplot, ecrplot, csrplot, fail, nocoord_old


# plot the load-deflection curve
def plot(fcVM, averaged, el_limit, ul_limit, un, lbd, csrplot, peeqmax, dl, du, target_LF, nstep, ue, ultimate_strain,
         disp_new, disp_old, elNodes, nocoord, sig_new, peeq, sigmises, csr, noce, sig_yield):
    class Index(object):

        def __init__(self, averaged, disp_new,
                     disp_old, elNodes, nocoord, sig_new, peeq,
                     sigmises, csr, noce, sig_yield):
            self.averaged = averaged
            self.elNodes = elNodes
            self.nocoord = nocoord
            self.sig_new = sig_new
            self.peeq = peeq
            self.sigmises = sigmises
            self.csr = csr
            self.noce = noce
            self.disp_old = disp_old
            self.disp_new = disp_new
            self.sig_yield = sig_yield
            self.ds = 1.0

        def stop(self, event):
            self.cnt = False
            plt.close()
            self.clicked = True

        def add(self, event):
            self.cnt = True
            self.clicked = True
            cr1 = (self.target_LF - self.LF) * (self.target_LF_out - self.LF) <= 0.0
            if (cr1):
                self.dl = np.sign(self.target_LF_out - self.LF) * 1.0 / self.nstep
                self.du = self.dl * self.ue
                print("self.dl: ", self.dl)
            plt.close()

        def rev(self, event):
            self.cnt = True
            self.dl = - self.dl
            self.du = - self.du
            plt.close()
            self.clicked = True

        def close_window(self, event):
            if self.cnt == False:
                self.stop('stop_event')

        def submit(self, LF):
            self.target_LF_out = float(LF)

        def set_scale(self, ds):
            self.ds = float(ds)

        def PSV(self, event):

            class Movie():
                def __init__(self, p):
                    pass

                def __call__(self, *args, **kwargs):
                    file = os.path.join(mdir, "output files", name + '_PSV.gif')
                    path = p.generate_orbital_path(n_points=36)
                    p.open_gif(file)
                    p.orbit_on_path(path, write_frames=True)

            class Toggle_Plane():
                def __init__(self, p):
                    pass

                def __call__(self, *args, **kwargs):
                    p.plane_widgets[0].SetEnabled(not p.plane_widgets[0].GetEnabled())

                def update(self):
                    pass

            class Scale_Stress():
                def __init__(self, p, scale):
                    self.scale = scale

                def __call__(self, value):
                    glyphs1 = grid1.glyph(orient="Major Principal Stress Vector", scale="Major Principal Stress",
                                          factor=10 ** value * self.scale,
                                          geom=geom)

                    p.add_mesh(glyphs1, name='sv1', show_scalar_bar=False, lighting=False, cmap=['red'])

                    glyphs2 = grid1.glyph(orient="Intermediate Principal Stress Vector",
                                          scale="Intermediate Principal Stress",
                                          factor=10 ** value * self.scale,
                                          geom=geom)

                    p.add_mesh(glyphs2, name='sv2', show_scalar_bar=False, lighting=False, cmap=['green'])

                    glyphs3 = grid1.glyph(orient="Minor Principal Stress Vector", scale="Minor Principal Stress",
                                          factor=10 ** value * self.scale,
                                          geom=geom)

                    p.add_mesh(glyphs3, name='sv3', show_scalar_bar=False, lighting=False, cmap=['blue'])

                def update(self):
                    pass

            def screen_shot(p):
                print(mdir)
                file = os.path.join(mdir, "output files", name + '_PSV.png')
                p.screenshot(file)

            # pv.global_theme.cmap = 'jet'
            pv.global_theme.cmap = 'coolwarm'
            # pv.global_theme.cmap = 'turbo'
            # pv.set_plot_theme('dark')
            pv.global_theme.full_screen = True
            pv.global_theme.title = 'PSV'

            # Controlling the text properties
            sargs = dict(
                title_font_size=16,
                label_font_size=14,
                shadow=False,
                color="white",
                bold=True,
                n_labels=5,
                italic=True,
                fmt="%.2e",
                font_family="arial",

            )

            tet10stress, tet10peeq, tet10csr, tet10svm, tet10triax = mapStresses(self.averaged, self.elNodes,
                                                                                 self.nocoord,
                                                                                 self.sig_new, self.peeq,
                                                                                 self.sigmises, self.csr, self.noce,
                                                                                 self.sig_yield)

            tet10s1, tet10s2, tet10s3, sv1, sv2, sv3 = calculate_principal_stress(tet10stress)

            x_range = max(nocoord[:, 0]) - min(nocoord[:, 0])
            y_range = max(nocoord[:, 1]) - min(nocoord[:, 1])
            z_range = max(nocoord[:, 2]) - min(nocoord[:, 2])

            geo_range = max(x_range, y_range, z_range)

            stress_range = max(max(map(abs, tet10s1)), max(map(abs, tet10s2)), max(map(abs, tet10s3)))

            if stress_range == 0.0:
                scale = 1.0
            else:
                scale = 0.2 * geo_range / stress_range

            disp_range = max(self.disp_new) - min(self.disp_new)

            # if disp_range == 0.0:
            #     scale_disp = 1.0
            # else:
            #     scale_disp = 0.2 * geo_range / disp_range

            padding = np.full(len(self.elNodes), 10, dtype=int)
            self.elements = np.vstack((padding, (self.elNodes - 1).T)).T

            celltypes = np.full(len(self.elNodes), fill_value=CellType.QUADRATIC_TETRA, dtype=np.uint8)
            geom = pv.Line()

            points = self.nocoord

            grid = pv.UnstructuredGrid(self.elements, celltypes, points)
            grid.point_data['Displacement'] = np.reshape(disp_new, (len(nocoord), 3))

            grid1 = grid.warp_by_vector(vectors='Displacement', factor=self.ds,
                                        inplace=False, progress_bar=False)
            grid1.point_data["Major Principal Stress"] = tet10s1.flatten(order="F")
            grid1.point_data["Intermediate Principal Stress"] = tet10s2.flatten(order="F")
            grid1.point_data["Minor Principal Stress"] = tet10s3.flatten(order="F")
            grid1.point_data['Major Principal Stress Vector'] = sv1
            grid1.point_data['Intermediate Principal Stress Vector'] = sv2
            grid1.point_data['Minor Principal Stress Vector'] = sv3

            glyphs1 = grid1.glyph(orient="Major Principal Stress Vector", scale="Major Principal Stress", factor=scale,
                                  geom=geom)
            glyphs2 = grid1.glyph(orient="Intermediate Principal Stress Vector", scale="Intermediate Principal Stress",
                                  factor=scale, geom=geom)
            glyphs3 = grid1.glyph(orient="Minor Principal Stress Vector", scale="Minor Principal Stress", factor=scale,
                                  geom=geom)

            p = pv.Plotter()

            p.set_background('grey', all_renderers=True)

            # p.set_background("royalblue", top="aliceblue")

            clip = p.add_mesh_clip_plane(grid1, name="mesh", show_edges=False, normal=[1.0, 0., 0.], invert=True,
                                         show_scalar_bar=False,
                                         scalar_bar_args=sargs, cmap=['grey'])
            p.add_mesh(glyphs1, name='sv1', show_scalar_bar=False, lighting=False, cmap=['red'])
            p.add_mesh(glyphs2, name='sv2', show_scalar_bar=False, lighting=False, cmap=['green'])
            p.add_mesh(glyphs3, name='sv3', show_scalar_bar=False, lighting=False, cmap=['blue'])

            p.add_key_event("s", lambda: screen_shot(p))

            ss = Scale_Stress(p, scale)
            p.add_slider_widget(
                callback=ss,
                rng=[-1.0, 1.0],
                value=0.0,
                title="log(Scale Stress)",
                pointa=(0.05, 0.075),
                pointb=(0.3, 0.075),
                style='modern',
                slider_width=0.02,
                tube_width=0.02
            )

            tp = Toggle_Plane(p)
            p.add_key_event("t", lambda: tp())

            p.show(cpos=[1.0, 1.0, 1.0])

        def VTK(self, event):

            class Toggle_Plane():
                def __init__(self, p):
                    self.normals = []
                    self.origins = []
                    self.names = ["mesh1", "mesh2", "mesh3", "mesh4"]

                def __call__(self, *args, **kwargs):
                    if p.plane_widgets:
                        self.normals = []
                        self.origins = []
                        for widget in p.plane_widgets:
                            self.normals.append(widget.GetNormal())
                            self.origins.append(widget.GetOrigin())
                        p.clear_plane_widgets()

                    else:
                        for name in self.names:
                            p.remove_actor(name)
                        p.subplot(0, 0)
                        p.add_mesh_clip_plane(grid1, name="mesh1", show_edges=False, normal=self.normals[0],
                                              origin=self.origins[0],
                                              invert=True, lighting=True,
                                              scalar_bar_args=sargs)
                        p.subplot(0, 1)
                        p.add_mesh_clip_plane(grid2, name="mesh2", show_edges=False, normal=self.normals[1],
                                              origin=self.origins[1],
                                              invert=True, lighting=True,
                                              scalar_bar_args=sargs)
                        p.subplot(1, 0)
                        p.add_mesh_clip_plane(grid3, name="mesh3", show_edges=False, normal=self.normals[2],
                                              origin=self.origins[2],
                                              invert=True, lighting=True,
                                              scalar_bar_args=sargs)
                        p.subplot(1, 1)
                        p.add_mesh_clip_plane(grid4, name="mesh4", show_edges=False, normal=self.normals[3],
                                              origin=self.origins[3],
                                              invert=True, lighting=True,
                                              scalar_bar_args=sargs)

                def update(self):
                    pass

            def screen_shot(p):
                print(mdir)
                file = os.path.join(mdir, "output files", name + '.png')
                p.screenshot(file)

            # pv.global_theme.cmap = 'jet'
            pv.global_theme.cmap = 'coolwarm'
            # pv.global_theme.cmap = 'turbo'
            # pv.set_plot_theme('dark')
            pv.global_theme.full_screen = True
            pv.global_theme.title = 'VTK'

            # Controlling the text properties
            sargs = dict(
                title_font_size=16,
                label_font_size=14,
                shadow=False,
                color="white",
                bold=True,
                n_labels=5,
                italic=True,
                fmt="%.2e",
                font_family="arial",

            )

            tet10stress, tet10peeq, tet10csr, tet10svm, tet10triax = mapStresses(self.averaged, self.elNodes,
                                                                                 self.nocoord,
                                                                                 self.sig_new, self.peeq,
                                                                                 self.sigmises, self.csr, self.noce,
                                                                                 self.sig_yield)

            x_range = max(nocoord[:, 0]) - min(nocoord[:, 0])
            y_range = max(nocoord[:, 1]) - min(nocoord[:, 1])
            z_range = max(nocoord[:, 2]) - min(nocoord[:, 2])

            geo_range = max(x_range, y_range, z_range)

            disp_range = max(self.disp_new) - min(self.disp_new)

            # if disp_range == 0.0:
            #     scale = 1.0
            # else:
            #     scale = self.ds * geo_range / disp_range

            padding = np.full(len(self.elNodes), 10, dtype=int)
            self.elements = np.vstack((padding, (self.elNodes - 1).T)).T

            celltypes = np.full(len(self.elNodes), fill_value=CellType.QUADRATIC_TETRA, dtype=np.uint8)

            points = self.nocoord + self.ds * np.reshape(disp_new, (len(nocoord), 3))
            grid1 = pv.UnstructuredGrid(self.elements, celltypes, points)
            grid1.point_data["Critical Strain Ratio\n"] = tet10csr.flatten(order="F")
            grid2 = pv.UnstructuredGrid(self.elements, celltypes, points)
            grid2.point_data["Equivalent Plastic Strain\n"] = tet10peeq.flatten(order="F")
            grid3 = pv.UnstructuredGrid(self.elements, celltypes, points)
            grid3.point_data["von Mises Stress\n"] = tet10svm.flatten(order="F")
            grid4 = pv.UnstructuredGrid(self.elements, celltypes, points)
            grid4.point_data["Triaxiality\n"] = tet10triax.flatten(order="F")

            p = pv.Plotter(shape=(2, 2))

            p.link_views()

            p.set_background('grey', all_renderers=True)

            # left upper pane
            p.subplot(0, 0)
            p.add_mesh_clip_plane(grid1, name="mesh1", show_edges=False, normal=[1.0, 0., 0.], invert=True,
                                  scalar_bar_args=sargs)

            # right upper pane
            p.subplot(0, 1)
            p.add_mesh_clip_plane(grid2, name="mesh2", show_edges=False, normal=[1.0, 0., 0.], invert=True,
                                  scalar_bar_args=sargs)

            # left lower pane
            p.subplot(1, 0)
            p.add_mesh_clip_plane(grid3, name="mesh3", show_edges=False, normal=[1.0, 0., 0.], invert=True,
                                  scalar_bar_args=sargs)

            # right lower pane
            p.subplot(1, 1)
            p.add_mesh_clip_plane(grid4, name="mesh4", show_edges=False, normal=[1.0, 0., 0.], invert=True,
                                  scalar_bar_args=sargs)

            p.add_key_event("s", lambda: screen_shot(p))

            tp = Toggle_Plane(p)
            p.add_key_event("t", lambda: tp())

            p.show(cpos=[1.0, 1.0, 1.0])

    print("displacement results: ", len(un), "load factor results: ", len(lbd))
    callback = Index(averaged, disp_new,
                     disp_old, elNodes, nocoord, sig_new, peeq,
                     sigmises, csr, noce, sig_yield)
    callback.cnt = False
    callback.clicked = False
    callback.dl = dl
    callback.du = du
    callback.target_LF = target_LF
    callback.target_LF_out = target_LF
    callback.LF = lbd[-1]
    callback.ue = ue
    callback.nstep = nstep
    fig, ax = plt.subplots(1, 2, figsize=(10, 6))
    ax[1].xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    fig.canvas.manager.set_window_title('fcVM')
    plt.subplots_adjust(bottom=0.2)
    ax[0].plot(un, lbd, '-ok')
    ax[0].set(xlabel='displacement [mm]', ylabel='load factor [-]', title='')
    if fcVM.csrRbtn.isChecked():
        ax[1].set(xlabel='critical strain ratio [-]', title='')
        ax[1].plot(csrplot, lbd, '-ok')
    else:
        ax[1].set(xlabel='equivalent plastic strain (PEEQ) [-]', title='')
        ax[1].plot(peeqmax, lbd, '-ok')
    ax[0].grid()
    ax[1].grid()
    b_w = 0.075
    b_h = 0.06
    b_s = 0.01
    b_y = 0.05
    axstop = plt.axes([0.3 - b_w / 2.0 - b_w - b_s, b_y, b_w, b_h])
    axadd = plt.axes([0.3 - b_w / 2.0, b_y, b_w, b_h])
    axbox = plt.axes([0.3 - b_w / 2.0 + b_w + b_s, b_y, b_w, b_h])
    dsbox = plt.axes([0.7 - b_w / 2.0 + b_w + b_s, b_y, b_w, b_h])
    axVTK = plt.axes([0.575, b_y, b_w, b_h])
    axPSV = plt.axes([0.575 + b_w + b_s, b_y, b_w, b_h])
    bstop = Button(axstop, 'stop')
    bstop.on_clicked(callback.stop)
    badd = Button(axadd, 'add')
    badd.on_clicked(callback.add)
    bVTK = Button(axVTK, 'VTK')
    bPSV = Button(axPSV, 'PSV')
    bVTK.on_clicked(callback.VTK)
    bPSV.on_clicked(callback.PSV)
    text_box = TextBox(axbox, "", textalignment="center")
    text_box.set_val(target_LF)
    text_box.on_submit(callback.submit)
    fig.text(0.3 - b_w / 2.0 + 2.0 * (b_w + b_s) - 0.005, b_y + b_h / 2.0 - 0.01, 'Target Load Factor', fontsize=10)
    fig.canvas.mpl_connect('close_event', callback.close_window)

    ds_box = TextBox(dsbox, "", textalignment="center")
    ds_box.set_val(1.0)
    ds_box.on_submit(callback.set_scale)
    fig.text(0.7 - b_w / 2.0 + 2.0 * (b_w + b_s) - 0.005, b_y + b_h / 2.0 - 0.01, 'Displacement Scale', fontsize=10)
    fig.canvas.mpl_connect('close_event', callback.close_window)


    if ul_limit != 0:
        if fcVM.csrRbtn.isChecked():
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
        else:
            fac = (ultimate_strain - peeqmax[ul_limit]) / (peeqmax[ul_limit + 1] - peeqmax[ul_limit])
            lbd_limit = lbd[ul_limit] + fac * (lbd[ul_limit + 1] - lbd[ul_limit])
            un_limit = un[ul_limit] + fac * (un[ul_limit + 1] - un[ul_limit])
            ax[0].plot([0.0, un_limit], [lbd_limit, lbd_limit], color='r', linestyle="--")
            ax[0].plot([un_limit, un_limit], [0.0, lbd_limit], color='r', linestyle="--")
            ax[0].plot([un[el_limit], un[el_limit]], [0.0, lbd[el_limit]], color='b', linestyle="--")
            ax[0].plot([0.0, un[el_limit]], [lbd[el_limit], lbd[el_limit]], color='b', linestyle="--")
            ax[1].plot([0.0, ultimate_strain], [lbd_limit, lbd_limit], color='r', linestyle="--")
            ax[1].plot([ultimate_strain, ultimate_strain], [0.0, lbd_limit], color='r', linestyle="--")
            ax[1].plot([0.0, ultimate_strain], [lbd[el_limit], lbd[el_limit]], color='b', linestyle="--")

    plt.show()

    while True:
        plt.pause(0.01)
        if callback.clicked:
            break

    return callback.cnt, callback.dl, callback.du, callback.target_LF_out


# update PEEQ and CSR
@jit(nopython=True, cache=True, fastmath=True)
def update_PEEQ_CSR(nelem, materialbyElement, sig_test, sig_new, sig_yield, ultimate_strain, peeq, csr, triax, pressure,
                    sigmises, ecr, Et_E):
    E = materialbyElement[0][0]  # Young's Modulus
    nu = materialbyElement[0][1]  # Poisson's Ratio
    G = E / 2.0 / (1 + nu)  # shear modulus
    if Et_E > 0.95: Et_E = 0.95
    Et = Et_E * E
    H = Et / (1.0 - Et_E)
    if ultimate_strain == 0.0:
        ultimate_strain = 1.0e12
    alpha = np.sqrt(np.e) * ultimate_strain  # stress triaxiality T = 1/3 for uniaxial test
    beta = 1.5

    for el in range(nelem):
        # print("sigma: ", sig_new[24 * el:24 * el + 6])
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

            DL = 0.0
            if (sig_mises_test > sig_yield[ipos1]):
                DL = (sig_mises_test - sig_yield[ipos1]) / (3.0 * G + H)
                peeq[ipos1] += DL
                sig_yield[ipos1] += Et * DL

            T = p_n / sig_yield[ipos1]  # experimental - verify against theory and tests

            pressure[ipos1] = p_n
            sigmises[ipos1] = sig_mises_new
            triax[ipos1] = T

            critical_strain = alpha * np.exp(-beta * T)

            if critical_strain < 1.0e-6:
                critical_strain = 1.0e-6

            ecr[ipos1] = critical_strain

            csr[ipos1] += DL / critical_strain


# convective stress update
@jit(nopython=True, cache=True, fastmath=True)
def convection(xlv, u, du, sig0, sigc, deps, dshp0, dshp1, dshp2, dmat):
    dshp = np.empty((10, 3), dtype=np.float64)  # xi, eta, zeta - derivative of the shape functions
    u10 = np.empty((10, 3), dtype=np.float64)  # xi, eta, zeta - derivative of the shape functions
    du10 = np.empty((10, 3), dtype=np.float64)  # xi, eta, zeta - derivative of the shape functions
    F = np.eye(3)

    for i in range(10):
        n3 = 3 * i
        u10[i] = u[n3:n3 + 3]
        du10[i] = du[n3:n3 + 3]
    dshp[:, 0] = dshp0
    dshp[:, 1] = dshp1
    dshp[:, 2] = dshp2

    sxx = sig0[0]
    syy = sig0[1]
    szz = sig0[2]
    sxy = sig0[3]
    szx = sig0[4]
    syz = sig0[5]

    sig = np.array([[sxx, sxy, szx],
                    [sxy, syy, syz],
                    [szx, syz, szz]])

    # local derivative of the global coordinates dx_j/dxi_m = dN_l/dxi_m * x0_jl (equation 7)
    xs = np.dot(xlv.T, dshp)

    # global derivative of the local coordinates dxi_m/dx_j = (dN_l/dxi_m * x0_jl)^(-1) (equation 7)
    xs_inv = np.linalg.inv(xs)

    # global derivatives of the shape functions (equation 6)
    dshpg = np.dot(dshp, xs_inv)

    # deformation gradient (equation 6)
    F = np.eye(3, 3, dtype=np.float64) + np.dot(du10.T, dshpg)  # deformation gradient

    # rho_1 / rho_0 (equation 3)
    rr = 1.0 / np.linalg.det(F)

    # stress update (equation 2)
    sig = rr * np.dot(F, np.dot(sig, F.T))

    sigc[0] = sig[0, 0]
    sigc[1] = sig[1, 1]
    sigc[2] = sig[2, 2]
    sigc[3] = sig[0, 1]
    sigc[4] = sig[0, 2]
    sigc[5] = sig[1, 2]

    return F


# update stresses and loads
@jit(nopython=True, cache=True, fastmath=True)
def update_stress_load(gp10, elNodes, nocoord, materialbyElement, sig_yield, disp_new, du, sig, sig_update,
                       sig_test_global, qin, Et_E, LD, pgp):
    du10x = np.empty(10, dtype=np.float64)  # displacement increments for the 10 tetrahedral nodes
    du10y = np.empty(10, dtype=np.float64)  # displacement increments for the 10 tetrahedral nodes
    du10z = np.empty(10, dtype=np.float64)  # displacement increments for the 10 tetrahedral nodes
    du10 = np.empty(30, dtype=np.float64)  # displacement increments for the 10 tetrahedral nodes
    u10 = np.empty(30, dtype=np.float64)  # displacement increments for the 10 tetrahedral nodes
    xlv = np.empty((10, 3), dtype=np.float64)
    sig_test = np.empty(6, dtype=np.float64)
    sigc = np.empty(6, dtype=np.float64)
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

    E = materialbyElement[0][0]  # Young's Modulus
    nu = materialbyElement[0][1]  # Poisson's Ratio
    G = E / 2.0 / (1 + nu)  # shear modulus
    if Et_E > 0.95: Et_E = 0.95
    Et = Et_E * E
    H = Et / (1.0 - Et_E)

    for el, nodes in enumerate(elNodes):
        elpos = 24 * el  # 4 integration points with 6 stress components each
        elv = np.zeros(30, dtype=np.float64)  # element load vector, 10 nodes with 3 load components each

        for index, nd in enumerate(nodes):
            n3 = 3 * (nd - 1)
            i3 = 3 * index
            du10x[index] = du[n3]
            du10y[index] = du[n3 + 1]
            du10z[index] = du[n3 + 2]
            du10[i3] = du[n3]
            du10[i3 + 1] = du[n3 + 1]
            du10[i3 + 2] = du[n3 + 2]
            u10[i3] = disp_new[n3]
            u10[i3 + 1] = disp_new[n3 + 1]
            u10[i3 + 2] = disp_new[n3 + 2]

        for index, nd in enumerate(nodes):
            i3 = 3 * index
            co = nocoord[nd - 1]
            if LD:
                xlv0[index] = co[0] + u10[i3]  # + du10[i3]
                xlv1[index] = co[1] + u10[i3 + 1]  # + du10[i3 + 1]
                xlv2[index] = co[2] + u10[i3 + 2]  # + du10[i3 + 2]
                xlv[index] = co + u10[i3:i3 + 3]
            else:
                xlv0[index] = co[0]
                xlv1[index] = co[1]
                xlv2[index] = co[2]
                xlv[index] = co

        for i in range(4):

            ipp = 4 * el + i
            sy = sig_yield[ipp]
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

            # strain
            deps = np.zeros(6, dtype=np.float64)
            for j in range(10):
                deps[0] += bm0[j] * du10x[j]
                deps[1] += bm1[j] * du10y[j]
                deps[2] += bm2[j] * du10z[j]
                deps[3] += bm3[j] * du10x[j] + bm4[j] * du10y[j]
                deps[4] += bm5[j] * du10x[j] + bm6[j] * du10z[j]
                deps[5] += bm7[j] * du10y[j] + bm8[j] * du10z[j]

            if LD:
                # F = convection(xlv, u10, du10, sig[ippos:ippos + 6], sigc, deps, dshp0, dshp1, dshp2, dmat)
                sxx = sig[ippos]
                syy = sig[ippos + 1]
                szz = sig[ippos + 2]
                sxy = sig[ippos + 3]
                szx = sig[ippos + 4]
                syz = sig[ippos + 5]

                st = np.array([[sxx, sxy, szx],
                               [sxy, syy, syz],
                               [szx, syz, szz]])

                F = np.eye(3, dtype=np.float64)  # deformation gradient
                for i in range(10):
                    F[0, 0] += du10x[i] * dshpg0[i]
                    F[0, 1] += du10x[i] * dshpg1[i]
                    F[0, 2] += du10x[i] * dshpg2[i]
                    F[1, 0] += du10y[i] * dshpg0[i]
                    F[1, 1] += du10y[i] * dshpg1[i]
                    F[1, 2] += du10y[i] * dshpg2[i]
                    F[2, 0] += du10z[i] * dshpg0[i]
                    F[2, 1] += du10z[i] * dshpg1[i]
                    F[2, 2] += du10z[i] * dshpg2[i]
                rr = (F[0, 0] * F[1, 1] * F[2, 2] -
                      F[0, 0] * F[1, 2] * F[2, 1] +
                      F[0, 2] * F[1, 0] * F[2, 1] -
                      F[0, 2] * F[1, 1] * F[2, 0] +
                      F[0, 1] * F[1, 2] * F[2, 0] -
                      F[0, 1] * F[1, 0] * F[2, 2])

                rr = 1.0 / rr

                sigcon = np.zeros((3, 3), dtype=np.float64)  # convected stress

                for i in range(3):
                    for j in range(3):
                        for k in range(3):
                            for l in range(3):
                                sigcon[i][k] += F[i][j] * st[j][l] * F[k][l]

                sigc[0] = rr * sigcon[0][0]
                sigc[1] = rr * sigcon[1][1]
                sigc[2] = rr * sigcon[2][2]
                sigc[3] = rr * sigcon[0][1]
                sigc[4] = rr * sigcon[0][2]
                sigc[5] = rr * sigcon[1][2]

            else:
                sigc = sig[ippos:ippos + 6]

            # elastic test stress
            for j in range(6):
                tmp = sigc[j]
                for k in range(6):
                    tmp += dmat[j, k] * deps[k]
                sig_test[j] = tmp

            sig_test_global[ippos:ippos + 6] = sig_test

            # return elastic stress to von Mises yield surface
            sxx, syy, szz, sxy, szx, syz, pp = vmises_original_optimised(sig_test, sy, H, G)
            sig_update[ippos:ippos + 6] = sxx, syy, szz, sxy, szx, syz
            pgp[ipp] = pp  # plastic Gauss point

            # calculate element load vectors
            ipxsj = ip[3] * abs(xsj)
            for j in range(10):
                j3 = 3 * j
                elv[j3] += (bm0[j] * sxx + bm3[j] * sxy + bm5[j] * szx) * ipxsj
                elv[j3 + 1] += (bm1[j] * syy + bm4[j] * sxy + bm7[j] * syz) * ipxsj
                elv[j3 + 2] += (bm2[j] * szz + bm6[j] * szx + bm8[j] * syz) * ipxsj

        # add element load vectors to the global load vector
        for i in range(10):
            iglob = nodes[i] - 1
            iglob3 = 3 * iglob
            i3 = 3 * i
            for k in range(3):
                qin[iglob3 + k] += elv[i3 + k]

    return


# @jit("types.UniTuple(float64, 6)(float64[::1], float64, float64, float64)", nopython=True, cache=True)
@jit(nopython=True, cache=True)
def vmises_original_optimised(sig_test, sig_yield, H, G):
    st0, st1, st2, st3, st4, st5 = sig_test
    p = (st0 + st1 + st2) / 3.0
    st0 -= p
    st1 -= p
    st2 -= p
    sig_mises = np.sqrt(1.5 * (st0 ** 2 + st1 ** 2 + st2 ** 2) +
                        3.0 * (st3 ** 2 + st4 ** 2 + st5 ** 2))

    if sig_yield > sig_mises:
        fac = 1.0
        pp = False  # elastic Gauss point
    else:
        fac = (1.0 - (1.0 - sig_yield / sig_mises) * 3.0 * G / (H + 3 * G))
        pp = True  # plastic Gauss point

    st0 = fac * st0 + p
    st1 = fac * st1 + p
    st2 = fac * st2 + p
    st3 = fac * st3
    st4 = fac * st4
    st5 = fac * st5

    return st0, st1, st2, st3, st4, st5, pp


# map stresses to nodes
@jit(nopython=True, cache=True)
def mapStresses(averaged, elNodes, nocoord, sig, peeq, sigvm, csr, noce, sig_yield):
    # map maps corner node stresses to all tet10 nodes

    map_inter = np.array([[0.5, 0.5, 0.0, 0.0],
                          [0.0, 0.5, 0.5, 0.0],
                          [0.5, 0.0, 0.5, 0.0],
                          [0.5, 0.0, 0.0, 0.5],
                          [0.0, 0.5, 0.0, 0.5],
                          [0.0, 0.0, 0.5, 0.5]])

    tet10stress = np.zeros((len(nocoord), 6), dtype=np.float64)
    tet10peeq = np.zeros(len(nocoord), dtype=np.float64)
    tet10csr = np.zeros(len(nocoord), dtype=np.float64)
    tet10svm = np.zeros(len(nocoord), dtype=np.float64)
    tet10triax = np.zeros(len(nocoord), dtype=np.float64)

    tet4stress = sig.reshape((len(elNodes), 4, 6))
    tet4peeq = peeq.reshape((len(elNodes), 4))
    tet4csr = csr.reshape((len(elNodes), 4))
    tet4svm = sigvm.reshape((len(elNodes), 4))
    tet4triax = (tet4stress[:, :, 0] + tet4stress[:, :, 1] + tet4stress[:, :, 2]) / 3.0 / sig_yield

    # averaged nodal stresses
    for el, nodes in enumerate(elNodes):
        corners = nodes[0:4]
        tet10stress[corners - 1] += tet4stress[el] / noce[corners - 1].reshape((4, 1))

    if averaged:
        # averaged nodal scalars
        for el, nodes in enumerate(elNodes):
            corners = nodes[0:4]
            tet10peeq[corners - 1] += tet4peeq[el] / noce[corners - 1]
            tet10csr[corners - 1] += tet4csr[el] / noce[corners - 1]
            tet10svm[corners - 1] += tet4svm[el] / noce[corners - 1]
            tet10triax[corners - 1] += tet4triax[el] / noce[corners - 1]
    else:
        # unaveraged nodal scalars
        for el, nodes in enumerate(elNodes):
            corners = nodes[0:4]
            tet10peeq[corners - 1] = np.fmax(tet10peeq[corners - 1], tet4peeq[el])
            tet10csr[corners - 1] = np.fmax(tet10csr[corners - 1], tet4csr[el])
            tet10svm[corners - 1] = np.fmax(tet10svm[corners - 1], tet4svm[el])
            tet10triax[corners - 1] = np.fmax(tet10triax[corners - 1], tet4triax[el])

    # for i in range(len(nocoord)):
    #     tet10triax[i] = (tet10stress[i][0] + tet10stress[i][1] + tet10stress[i][2]) / 3.0 / sig_yield

    # results intermediate nodes
    for el, nodes in enumerate(elNodes):
        nd_corner = nodes[0:4]
        nd_inter = nodes[4:10]
        tet10stress[nd_inter - 1] = np.dot(map_inter, tet10stress[nd_corner - 1])
        tet10peeq[nd_inter - 1] = np.dot(map_inter, tet10peeq[nd_corner - 1])
        tet10csr[nd_inter - 1] = np.dot(map_inter, tet10csr[nd_corner - 1])
        tet10svm[nd_inter - 1] = np.dot(map_inter, tet10svm[nd_corner - 1])
        tet10triax[nd_inter - 1] = np.dot(map_inter, tet10triax[nd_corner - 1])

    return tet10stress, tet10peeq, tet10csr, tet10svm, tet10triax


# fill resultobject with results
def pasteResults(gnl, doc, elNodes, nocoord, nocoord_old, dis, disp_el, eigenval, eigenvec, tet10stress, tet10peeq,
                 tet10csr, tet10svm, maxImp, nstep):
    maxImp = float(maxImp)
    return_code = 0
    for obj in doc.Objects:
        if obj.Name[:8] == "Analysis":
            analysis = doc.getObject(obj.Name)

    if gnl == "GNLY" and not (float(nstep) > 1.0 and maxImp == 0.0):
        evec1 = eigenvec[:, 0] / np.max(np.abs(eigenvec[:, 0]))
        evec2 = eigenvec[:, 1] / np.max(np.abs(eigenvec[:, 1]))

    # analysis = doc.getObject("Analysis")

    nn = len(nocoord)  # number of nodes

    elements_tetra10 = {}
    mode_results_vol = {}
    mode_results_vol2 = {}
    mode_results_vol3 = {}
    mode_results_vol4 = {}

    nodes = range(1, nn + 1)
    coordinates = map(App.Vector, nocoord)
    coordinates_old = map(App.Vector, nocoord_old)
    volnodes = dict(zip(nodes, coordinates))
    volnodes_old = dict(zip(nodes, coordinates_old))

    resVol = analysis.addObject(ObjectsFem.makeResultMechanical(doc))[0]
    resVol.Label = "ElasticPlasticDisplacementResults"
    displacements = map(App.Vector, np.array_split(dis, nn))
    mode_disp_vol = dict(zip(nodes, displacements))
    mode_results_vol['disp'] = mode_disp_vol
    results = [mode_results_vol]

    resVol2 = analysis.addObject(ObjectsFem.makeResultMechanical(doc))[0]
    resVol2.Label = "ElasticDisplacementResults"
    displacements2 = map(App.Vector, np.array_split(disp_el, nn))
    mode_disp_vol2 = dict(zip(nodes, displacements2))
    mode_results_vol2['disp'] = mode_disp_vol2
    results2 = [mode_results_vol2]

    if gnl == "GNLY" and not (float(nstep) > 1.0 and maxImp == 0.0):
        resVol3 = analysis.addObject(ObjectsFem.makeResultMechanical(doc))[0]
        resVol3.Label = "ElasticBucklingShape_lambda1=" + str(round(eigenval[0], 3))
        displacements3 = map(App.Vector, np.array_split(evec1, nn))
        mode_disp_vol3 = dict(zip(nodes, displacements3))
        mode_results_vol3['disp'] = mode_disp_vol3
        results3 = [mode_results_vol2]

        resVol4 = analysis.addObject(ObjectsFem.makeResultMechanical(doc))[0]
        resVol4.Label = "ElasticBucklingShape_lambda2=" + str(round(eigenval[1], 3))
        displacements4 = map(App.Vector, np.array_split(evec2, nn))
        mode_disp_vol4 = dict(zip(nodes, displacements4))
        mode_results_vol4['disp'] = mode_disp_vol4
        results4 = [mode_results_vol2]

    for index, elem in enumerate(elNodes):
        elements_tetra10[index + 1] = (
            elem[0], elem[2], elem[1], elem[3], elem[6], elem[5], elem[4], elem[7], elem[9], elem[8])

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

    mvol2 = {
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
        'Results': results2
    }

    if gnl == "GNLY" and not (float(nstep) > 1.0 and maxImp == 0.0):
        mvol3 = {
            'Nodes': volnodes_old,
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
            'Results': results3
        }

        mvol4 = {
            'Nodes': volnodes_old,
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
            'Results': results4
        }

    meshvol = itf.make_femmesh(mvol)
    meshvol2 = itf.make_femmesh(mvol2)
    if gnl == "GNLY" and not (float(nstep) > 1.0 and maxImp == 0.0):
        meshvol3 = itf.make_femmesh(mvol3)
        meshvol4 = itf.make_femmesh(mvol4)

    result_mesh_object_1 = ObjectsFem.makeMeshResult(doc, 'Result_Mesh_Volume')
    result_mesh_object_1.FemMesh = meshvol

    result_mesh_object_2 = ObjectsFem.makeMeshResult(doc, 'Result_Mesh_Volume')
    result_mesh_object_2.FemMesh = meshvol2

    if gnl == "GNLY" and not (float(nstep) > 1.0 and maxImp == 0.0):
        result_mesh_object_3 = ObjectsFem.makeMeshResult(doc, 'Result_Mesh_Volume')
        result_mesh_object_3.FemMesh = meshvol3

        result_mesh_object_4 = ObjectsFem.makeMeshResult(doc, 'Result_Mesh_Volume')
        result_mesh_object_4.FemMesh = meshvol4

    resVol.DisplacementVectors = [App.Vector(dis[3 * n], dis[3 * n + 1], dis[3 * n + 2]) for n in range(nn)]
    resVol.DisplacementLengths = [np.linalg.norm([dis[3 * n], dis[3 * n + 1], dis[3 * n + 2]]) for n in range(nn)]
    resVol.NodeStressXX = tet10stress.T[0].T.tolist()
    resVol.NodeStressYY = tet10stress.T[1].T.tolist()
    resVol.NodeStressZZ = tet10stress.T[2].T.tolist()
    resVol.NodeStressXY = tet10stress.T[3].T.tolist()
    resVol.NodeStressXZ = tet10stress.T[4].T.tolist()
    resVol.NodeStressYZ = tet10stress.T[5].T.tolist()

    resVol2.DisplacementVectors = [App.Vector(disp_el[3 * n], disp_el[3 * n + 1], disp_el[3 * n + 2]) for n in
                                   range(nn)]
    resVol2.DisplacementLengths = [np.linalg.norm([disp_el[3 * n], disp_el[3 * n + 1], disp_el[3 * n + 2]]) for n in
                                   range(nn)]

    if gnl == "GNLY" and not (float(nstep) > 1.0 and maxImp == 0.0):
        resVol3.DisplacementVectors = [
            App.Vector(evec1[3 * n], evec1[3 * n + 1], evec1[3 * n + 2]) for n in
            range(nn)]
        resVol3.DisplacementLengths = [
            np.linalg.norm([evec1[3 * n], evec1[3 * n + 1], evec1[3 * n + 2]]) for n in
            range(nn)]
        resVol4.DisplacementVectors = [
            App.Vector(evec2[3 * n], evec2[3 * n + 1], evec2[3 * n + 2]) for n in range(nn)]
        resVol4.DisplacementLengths = [
            np.linalg.norm([evec2[3 * n], evec2[3 * n + 1], evec2[3 * n + 2]]) for n in
            range(nn)]

    try:
        resVol.CriticalStrainRatio = tet10csr.T.tolist()  # FreeCAD 0.21.0 and higher - works for export to VTK
        resVol.Peeq = tet10peeq.T.tolist()
        resVol.Temperature = tet10csr.T.tolist()
        resVol.vonMises = tet10svm.T.tolist()

    except:
        resVol.Peeq = tet10csr.T.tolist()  # a hack for FreeCAD 0.20.x - store the critical strain ratio in the PEEQ output
        resVol.Temperature = tet10csr.T.tolist()
        resVol.vonMises = tet10svm.T.tolist()

    resVol.Stats = [min(dis[0::3]), max(dis[0::3]),
                    min(dis[1::3]), max(dis[1::3]),
                    min(dis[2::3]), max(dis[2::3]),
                    min(resVol.DisplacementLengths), max(resVol.DisplacementLengths),
                    min(resVol.vonMises), max(resVol.vonMises),
                    0.0, 0.0,
                    0.0, 0.0,
                    0.0, 0.0,
                    0.0, 0.0,
                    min(resVol.Peeq), max(resVol.Peeq),
                    min(resVol.Temperature), max(resVol.Temperature),
                    0.0, 0.0,
                    0.0, 0.0]

    resVol2.Stats = [min(disp_el[0::3]), max(disp_el[0::3]),
                     min(disp_el[1::3]), max(disp_el[1::3]),
                     min(disp_el[2::3]), max(disp_el[2::3]),
                     min(resVol2.DisplacementLengths), max(resVol2.DisplacementLengths),
                     0.0, 0.0,
                     0.0, 0.0,
                     0.0, 0.0,
                     0.0, 0.0,
                     0.0, 0.0,
                     0.0, 0.0,
                     0.0, 0.0,
                     0.0, 0.0,
                     0.0, 0.0]

    if gnl == "GNLY" and not (float(nstep) > 1.0 and maxImp == 0.0):
        resVol3.Stats = [min(evec1[0::3]), max(evec1[0::3]),
                         min(evec1[1::3]), max(evec1[1::3]),
                         min(evec1[2::3]), max(evec1[2::3]),
                         min(resVol3.DisplacementLengths), max(resVol3.DisplacementLengths),
                         0.0, 0.0,
                         0.0, 0.0,
                         0.0, 0.0,
                         0.0, 0.0,
                         0.0, 0.0,
                         0.0, 0.0,
                         0.0, 0.0,
                         0.0, 0.0,
                         0.0, 0.0]
        resVol4.Stats = [min(evec2[0::3]), max(evec2[0::3]),
                         min(evec2[1::3]), max(evec2[1::3]),
                         min(evec2[2::3]), max(evec2[2::3]),
                         min(resVol4.DisplacementLengths), max(resVol4.DisplacementLengths),
                         0.0, 0.0,
                         0.0, 0.0,
                         0.0, 0.0,
                         0.0, 0.0,
                         0.0, 0.0,
                         0.0, 0.0,
                         0.0, 0.0,
                         0.0, 0.0,
                         0.0, 0.0]

    resVol.Mesh = result_mesh_object_1
    resVol.NodeNumbers = [int(key) for key in resVol.Mesh.FemMesh.Nodes.keys()]

    resVol2.Mesh = result_mesh_object_2
    resVol2.NodeNumbers = [int(key) for key in resVol2.Mesh.FemMesh.Nodes.keys()]

    if gnl == "GNLY" and not (float(nstep) > 1.0 and maxImp == 0.0):
        resVol3.Mesh = result_mesh_object_3
        resVol3.NodeNumbers = [int(key) for key in resVol3.Mesh.FemMesh.Nodes.keys()]

        resVol4.Mesh = result_mesh_object_4
        resVol4.NodeNumbers = [int(key) for key in resVol4.Mesh.FemMesh.Nodes.keys()]

    resVol = itf.fill_femresult_mechanical(resVol, results)

    itf.fill_femresult_mechanical(resVol2, results2)

    if gnl == "GNLY" and not (float(nstep) > 1.0 and maxImp == 0.0):
        itf.fill_femresult_mechanical(resVol3, results3)
        itf.fill_femresult_mechanical(resVol4, results4)

    trm._TaskPanel.result_obj = resVol
    trm._TaskPanel.mesh_obj = meshvol

    doc.recompute()

    return resVol


def calcSum(Edge_Elements, Face_Elements, mesh, CSR, peeq, svm):
    gp10, gp6, gp2 = gaussPoints()

    coor = mesh.Nodes

    edge_length = []
    edge_peeq = []
    edge_CSR = []
    edge_svm = []

    for index, edge in enumerate(Edge_Elements):
        edge_length.append(0.0)
        edge_peeq.append(0.0)
        edge_CSR.append(0.0)
        edge_svm.append(0.0)
        for element in edge:
            xle = []
            for node in element:
                xle.append([coor[node].x, coor[node].y, coor[node].z])
            # integrate variable
            xle = np.array(xle).T
            for gp in gp2:
                xi = gp[0]
                xsj, shp = shape2lin(xi, xle)
                for i in range(3):
                    nd = element[i] - 1
                    dl = shp[i] * abs(xsj) * gp[1]
                    edge_length[index] += dl
                    edge_peeq[index] += peeq[nd] * dl
                    edge_CSR[index] += CSR[nd] * dl
                    edge_svm[index] += svm[nd] * dl
        Length = edge_length[index]
        if Length > 0.0:
            edge_peeq[index] /= Length
            edge_CSR[index] /= Length
            edge_svm[index] /= Length

    face_area = []
    face_peeq = []
    face_CSR = []
    face_svm = []

    for index, face in enumerate(Face_Elements):
        face_area.append(0.0)
        face_peeq.append(0.0)
        face_CSR.append(0.0)
        face_svm.append(0.0)
        for element in face:
            xlf = []
            for node in element:
                xlf.append([coor[node].x, coor[node].y, coor[node].z])
            # integrate variable
            xlf = np.array(xlf).T
            for gp in gp6:
                xi = gp[0]
                et = gp[1]
                xsj, shp, bmatS, xx, xt, xp = shape6tri(xi, et, xlf)
                for i in range(6):
                    nd = element[i] - 1
                    dA = shp[i] * abs(xsj) * gp[2]
                    face_area[index] += dA
                    face_peeq[index] += peeq[nd] * dA
                    face_CSR[index] += CSR[nd] * dA
                    face_svm[index] += svm[nd] * dA
        Area = face_area[index]
        if Area > 0.0:
            face_peeq[index] /= Area
            face_CSR[index] /= Area
            face_svm[index] /= Area

    return edge_length, edge_peeq, edge_CSR, edge_svm, face_area, face_peeq, face_CSR, face_svm


def exportVTK(gnl, elNodes, nocoord, nocoord_old, dis, disp_el, eigenval, eigenvec, tet10stress, tet10peeq, tet10csr,
              tet10svm, tet10triax, fy, file, maxImp, nstep):
    maxImp = float(maxImp)

    padding = np.full(len(elNodes), 10, dtype=int)
    cells = np.vstack((padding, (elNodes - 1).T)).T

    celltypes = np.full(len(elNodes), fill_value=CellType.QUADRATIC_TETRA, dtype=np.uint8)

    grid = pv.UnstructuredGrid(cells, celltypes, nocoord)
    grid_old = pv.UnstructuredGrid(cells, celltypes, nocoord)
    grid.point_data["Critical Strain Ratio\n"] = tet10csr.flatten(order="F")
    grid.point_data["Equivalent Plastic Strain\n"] = tet10peeq.flatten(order="F")
    grid.point_data["von Mises Stress\n"] = tet10svm.flatten(order="F")
    grid.point_data["Triaxiality\n"] = tet10triax.flatten(order="F")

    displacement = dis.reshape((len(nocoord), 3))
    elastic_disp = disp_el.reshape((len(nocoord), 3))
    stress = tet10stress.reshape((len(nocoord), 6))

    grid.point_data['Displacement'] = displacement
    if gnl == "GNLY" and not (float(nstep) > 1.0 and maxImp == 0.0):
        evec1 = eigenvec[:, 0] / np.max(np.abs(eigenvec[:, 0]))
        evec2 = eigenvec[:, 1] / np.max(np.abs(eigenvec[:, 1]))
        evec1 = evec1.reshape((len(nocoord), 3))
        evec2 = evec2.reshape((len(nocoord), 3))
        grid.point_data['Elastic Displacement'] = elastic_disp
        grid.point_data['Buckling shape for lambda1 = ' + str(round(eigenval[0], 3))] = evec1
        grid.point_data['Buckling shape for lambda2 = ' + str(round(eigenval[1], 3))] = evec2

    grid.point_data['Stress Tensor'] = stress

    tet10s1, tet10s2, tet10s3, sv1, sv2, sv3 = calculate_principal_stress(tet10stress)

    tet10rho = calculate_rho(tet10stress, fy)

    grid.point_data["Major Principal Stress\n"] = tet10s1.flatten(order="F")
    grid.point_data["Intermediate Principal Stress\n"] = tet10s2.flatten(order="F")
    grid.point_data["Minor Principal Stress\n"] = tet10s3.flatten(order="F")
    grid.point_data['Major Principal Stress Vector'] = sv1
    grid.point_data['Intermediate Principal Stress Vector'] = sv2
    grid.point_data['Minor Principal Stress Vector'] = sv3

    # grid.point_data['Reinforcement Ratio x'] = tet10rho[:, 0].flatten(order="F")
    # grid.point_data['Reinforcement Ratio y'] = tet10rho[:, 1].flatten(order="F")
    # grid.point_data['Reinforcement Ratio z'] = tet10rho[:, 2].flatten(order="F")

    pv.save_meshio(file, grid)


@jit(nopython=True, cache=True, nogil=True)
def calculate_principal_stress(tet10stress):
    s1 = np.zeros((len(tet10stress), 3), dtype=np.float64)
    s2 = np.zeros((len(tet10stress), 3), dtype=np.float64)
    s3 = np.zeros((len(tet10stress), 3), dtype=np.float64)
    tet10s1 = np.zeros(len(tet10stress), dtype=np.float64)
    tet10s2 = np.zeros(len(tet10stress), dtype=np.float64)
    tet10s3 = np.zeros(len(tet10stress), dtype=np.float64)

    for index, stress in enumerate(tet10stress):
        s11 = stress[0]  # Sxx
        s22 = stress[1]  # Syy
        s33 = stress[2]  # Szz
        s12 = stress[3]  # Sxy
        s31 = stress[4]  # Szx
        s23 = stress[5]  # Syz
        sigma = np.array([
            [s11, s12, s31],
            [s12, s22, s23],
            [s31, s23, s33]
        ])

        # print(sigma)

        eigenvalues, eigenvectors = np.linalg.eigh(sigma)

        eigenvalues = eigenvalues.real
        eigenvectors = eigenvectors.real

        idx = eigenvalues.argsort()[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]

        tet10s1[index] = eigenvalues[0]
        tet10s2[index] = eigenvalues[1]
        tet10s3[index] = eigenvalues[2]

        s1[index] = eigenvalues[0] * eigenvectors[:, 0]
        s2[index] = eigenvalues[1] * eigenvectors[:, 1]
        s3[index] = eigenvalues[2] * eigenvectors[:, 2]

    return tet10s1, tet10s2, tet10s3, s1, s2, s3


def calculate_rho(tet10stress, fy):
    #
    #   Calculation of Reinforcement Ratios and
    #   Concrete Stresses according to http://heronjournal.nl/53-4/3.pdf
    #           - See post:
    #             https://forum.freecadweb.org/viewtopic.php?f=18&t=28821
    #                   fy: factored yield strength of reinforcement bars
    #

    tet10rho = np.zeros((len(tet10stress), 3), dtype=np.float64)

    for index, stress_tensor in enumerate(tet10stress):
        rmin = 1.0e9
        eqmin = 14

        sxx = stress_tensor[0]
        syy = stress_tensor[1]
        szz = stress_tensor[2]
        sxy = stress_tensor[3]
        syz = stress_tensor[5]
        sxz = stress_tensor[4]

        rhox = np.zeros(15)
        rhoy = np.zeros(15)
        rhoz = np.zeros(15)

        #    i1=sxx+syy+szz NOT USED
        #    i2=sxx*syy+syy*szz+szz*sxx-sxy**2-sxz**2-syz**2 NOT USED
        i3 = (sxx * syy * szz + 2 * sxy * sxz * syz - sxx * syz ** 2
              - syy * sxz ** 2 - szz * sxy ** 2)

        #    Solution (5)
        d = (sxx * syy - sxy ** 2)
        if d != 0.:
            rhoz[0] = i3 / d / fy

        #    Solution (6)
        d = (sxx * szz - sxz ** 2)
        if d != 0.:
            rhoy[1] = i3 / d / fy

        #    Solution (7)
        d = (syy * szz - syz ** 2)
        if d != 0.:
            rhox[2] = i3 / d / fy

        #    Solution (9)
        if sxx != 0.:
            fc = sxz * sxy / sxx - syz
            fxy = sxy ** 2 / sxx
            fxz = sxz ** 2 / sxx

            #    Solution (9+)
            rhoy[3] = syy - fxy + fc
            rhoy[3] /= fy
            rhoz[3] = szz - fxz + fc
            rhoz[3] /= fy

            #    Solution (9-)
            rhoy[4] = syy - fxy - fc
            rhoy[4] /= fy
            rhoz[4] = szz - fxz - fc
            rhoz[4] /= fy

        #   Solution (10)
        if syy != 0.:
            fc = syz * sxy / syy - sxz
            fxy = sxy ** 2 / syy
            fyz = syz ** 2 / syy

            # Solution (10+)
            rhox[5] = sxx - fxy + fc
            rhox[5] /= fy
            rhoz[5] = szz - fyz + fc
            rhoz[5] /= fy

            # Solution (10-)vm
            rhox[6] = sxx - fxy - fc

            rhox[6] /= fy
            rhoz[6] = szz - fyz - fc
            rhoz[6] /= fy

        # Solution (11)
        if szz != 0.:
            fc = sxz * syz / szz - sxy
            fxz = sxz ** 2 / szz
            fyz = syz ** 2 / szz

            # Solution (11+)
            rhox[7] = sxx - fxz + fc
            rhox[7] /= fy
            rhoy[7] = syy - fyz + fc
            rhoy[7] /= fy

            # Solution (11-)
            rhox[8] = sxx - fxz - fc
            rhox[8] /= fy
            rhoy[8] = syy - fyz - fc
            rhoy[8] /= fy

        # Solution (13)
        rhox[9] = (sxx + sxy + sxz) / fy
        rhoy[9] = (syy + sxy + syz) / fy
        rhoz[9] = (szz + sxz + syz) / fy

        # Solution (14)
        rhox[10] = (sxx + sxy - sxz) / fy
        rhoy[10] = (syy + sxy - syz) / fy
        rhoz[10] = (szz - sxz - syz) / fy

        # Solution (15)
        rhox[11] = (sxx - sxy - sxz) / fy
        rhoy[11] = (syy - sxy + syz) / fy
        rhoz[11] = (szz - sxz + syz) / fy

        # Solution (16)
        rhox[12] = (sxx - sxy + sxz) / fy
        rhoy[12] = (syy - sxy - syz) / fy
        rhoz[12] = (szz + sxz - syz) / fy

        # Solution (17)
        if syz != 0.:
            rhox[13] = (sxx - sxy * sxz / syz) / fy
        if sxz != 0.:
            rhoy[13] = (syy - sxy * syz / sxz) / fy
        if sxy != 0.:
            rhoz[13] = (szz - sxz * syz / sxy) / fy

        for ir in range(0, rhox.size):

            if rhox[ir] >= -1.e-10 and rhoy[ir] >= -1.e-10 and rhoz[ir] > -1.e-10:

                # Concrete Stresses
                scxx = sxx - rhox[ir] * fy
                scyy = syy - rhoy[ir] * fy
                sczz = szz - rhoz[ir] * fy
                ic1 = (scxx + scyy + sczz)
                ic2 = (scxx * scyy + scyy * sczz + sczz * scxx - sxy ** 2
                       - sxz ** 2 - syz ** 2)
                ic3 = (scxx * scyy * sczz + 2 * sxy * sxz * syz - scxx * syz ** 2
                       - scyy * sxz ** 2 - sczz * sxy ** 2)

                if ic1 <= 1.e-6 and ic2 >= -1.e-6 and ic3 <= 1.0e-6:

                    rsum = rhox[ir] + rhoy[ir] + rhoz[ir]

                    if rsum < rmin and rsum > 0.:
                        rmin = rsum
                        eqmin = ir

        tet10rho[index] = [rhox[eqmin], rhoy[eqmin], rhoz[eqmin]]

    return tet10rho


def calculate_mohr_coulomb(prin1, prin3, phi, fck):
    #
    #             Calculation of Mohr Coulomb yield criterion to judge
    #             concrete curshing and shear failure
    #                   phi: angle of internal friction
    #                   fck: factored compressive strength of the matrix material (usually concrete)
    #

    coh = fck * (1 - np.sin(phi)) / 2 / np.cos(phi)

    mc_stress = ((prin1 - prin3) + (prin1 + prin3) * np.sin(phi)
                 - 2. * coh * np.cos(phi))

    if mc_stress < 0.:
        mc_stress = 0.

    return mc_stress
