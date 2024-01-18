# ***************************************************************************
# *                                                                         *
# *   Copyright (c) 2019 - Harry van Langen <hvlanalysis@icloud.com>        *
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


import Fem
import sys
import math
import time
import ObjectsFem
import numpy as np
import Part as Part
import FreeCAD as App
import FreeCADGui as Gui
from FreeCAD import Units
import scipy.sparse as scsp
from numba.typed import Dict
import matplotlib.pyplot as plt
from femtools import membertools
from femmesh import meshsetsgetter
from femmesh import meshtools as mt
from numba import types, __version__
from sksparse.cholmod import cholesky
from matplotlib.widgets import Button
from numba import jit, prange, objmode
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
    femmesh = docmesh.FemMesh
    member = membertools.AnalysisMember(analysis)

    # determine elements connected to a node using FC API
    fet = mt.get_femelement_table(mesh)
    # fet is dictionary: { elementid : [ nodeid, nodeid, ... , nodeid ] }
    net = mt.get_femnodes_ele_table(mesh.Nodes, fet)
    # net is dictionary: {nodeID : [[eleID, binary node position], [], ...], nodeID : [[], [], ...], ...}
    # node0 has binary node position 2^0 = 1, node1 = 2^1 = 2, ..., node10 = 2^10 = 1024

    # create connectivity array elNodes for mapping local node number -> global node number
    elNodes = np.array([mesh.getElementNodes(el) for el in mesh.Volumes])  # elNodes[elementIndex] = [node1,...,Node10]
    elo = dict(zip(mesh.Volumes, range(len(mesh.Volumes))))  # elo : {elementNumber : elementIndex}

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
    ppEl = {}  # BooleanFragment Primitive element El belonga to
    materialbyElement = []  # see further

    prn_upd("Number of material objects: ", len(member.mats_linear))

    for obj in App.ActiveDocument.Objects:
        if obj.Name == "BooleanFragments":
            num_BF_els = len(App.ActiveDocument.BooleanFragments.Objects)  # number of primitives in BooleanFragments
            num_Mat_obs = len(member.mats_linear)  # number of material objects
            if num_BF_els != num_Mat_obs:
                prn_upd("Each BooleanFragment primitive needs its own material object")
                raise SystemExit()
            for indp, primobject in enumerate(obj.Objects):  # primobject: primitive in BooleanFragment
                for indm, matobject in enumerate(member.mats_linear):  # matobject: material object
                    if matobject['Object'].References[0][0] == primobject:
                        matCon[indm] = indp  # the BooleanFragment primitive the material object refers to

    for indm, matobject in enumerate(member.mats_linear):
        E = float(App.Units.Quantity(matobject['Object'].Material['YoungsModulus']).getValueAs('MPa'))
        Nu = float(matobject['Object'].Material['PoissonRatio'])
        prn_upd("Material Object: ", matobject['Object'].Name, "   E= ", E, "   Nu= ", Nu)
        for el in element_sets[indm]:  # element_sets[indm]: all elements with material indm
            if matCon: ppEl[el] = matCon[indm]  # ppEl[el]: primitive el belongs to
            materialbyElement.append([E, Nu])  # materialbyElement[elementIndex] = [E, Nu]

    materialbyElement = np.asarray(materialbyElement)

    # set up interface element connectivity
    nodecount = len(nocoord)
    interface_elements = []
    tk = []
    tv = []
    ti = []
    NiP = {}
    shiftelements = []

    for obj in App.ActiveDocument.Objects:
        if obj.Name == "BooleanFragments":
            for indp1, primitive1 in enumerate(obj.Objects):  # primitive1: n-th primitive in BooleanFragment
                for face1 in primitive1.Shape.Faces:  # face: a face of n-th primitive
                    elfaces1 = np.asarray(
                        mesh.getFacesByFace(face1))  # el1faces[primitive face index] = [6-node-face numbers]
                    for sixNdFace in elfaces1:  # sixNdFace: 6-node face element
                        nodes = mesh.getElementNodes(sixNdFace)
                        for nd in nodes:
                            if nd not in NiP.keys():  # NiP[node] = Primitive Number
                                NiP[nd] = indp1
                    for indp2, primitive2 in enumerate(obj.Objects[indp1 + 1:]):
                        for face2 in primitive2.Shape.Faces:  # face: a face of n+1-th primitive
                            elfaces2 = np.asarray(
                                mesh.getFacesByFace(face2))  # el1faces[primitive face index] = [6-node-face numbers]
                            contact = list(set(elfaces1) & set(elfaces2))  # triangles in contact
                            if contact:
                                npflist = []  # contact nodes on primitive face
                                for sixNdFace in elfaces1:  # sixNdFace: 6-node-face element in face
                                    nodes = mesh.getElementNodes(sixNdFace)
                                    npflist.append(nodes)  # add nodes of sixNdFace to primitive contact face
                                cn = list(set([node for npf in npflist for node in
                                               npf]))  # remove duplicates
                                cn.sort()  # sort nodes
                                cn1_new = list(range(nodecount + 1, nodecount + len(cn) + 1))  # new nodes on face n
                                cn2_new = list(range(nodecount + len(cn) + 1,
                                                     nodecount + 2 * len(cn) + 1))  # new nodes on face n+1
                                nodecount += 2 * len(cn)
                                tk += cn + cn  # reference node numbers
                                tv += cn1_new + cn2_new  # new node numbers
                                ti += len(cn) * [indp1] + len(cn) * [
                                    indp1 + indp2 + 1]  # parent primitive number of new node
                                for i in range(2):  # for 2 sets of new nodes
                                    for node in cn:
                                        nocoord = np.append(nocoord, [nocoord[node - 1]],
                                                            axis=0)  # add coordinates for new nodes
                                for sixNdFace in contact:
                                    rn = mesh.getElementNodes(sixNdFace)  # reference nodes
                                    nodeset1 = []
                                    nodeset2 = []
                                    for nd in rn:
                                        nodeset1.append(cn1_new[cn.index(nd)])
                                        nodeset2.append(cn2_new[cn.index(nd)])
                                    interface_elements.append(
                                        nodeset1 + nodeset2)  # interface_elements[index] = [nodeA1,..,nodeA6, nodeB1,..,nodeB6]
                                for node in cn:
                                    for el in net[
                                        node]:  # net: {nodeID : [[eleID, binary node position], [], ...], nodeID : ...}
                                        nodepos = el[1]  # binary node position
                                        nonum = int(math.log(nodepos,
                                                             2))  # local node number from binary node number net[node][1]
                                        if ppEl[el[0]] == indp1:  # if element embedded in primitive n
                                            new_node = cn1_new[cn.index(node)]
                                            elNodes[elo[el[0]]][
                                                nonum] = new_node  # connect element to new node in primitive n
                                            NiP[new_node] = indp1
                                        elif ppEl[el[0]] == indp1 + indp2 + 1:  # if element embedded in primitive n+1
                                            new_node = cn2_new[cn.index(node)]
                                            elNodes[elo[el[0]]][
                                                nonum] = new_node  # connect element to new node in primitive n+1
                                            NiP[new_node] = indp1 + indp2 + 1
                                        shiftelements.append(elo[el[0]])

    twins = {}  # twins: {ref node 1: {primitive 1: [new_node 1, new_node 2], etc.}, etc.}
    for ind, node in enumerate(tk):  # for reference nodes
        if node not in twins:
            twins[node] = {}
        if ti[ind] not in twins[node]:
            twins[node][ti[ind]] = []
        twins[node][ti[ind]].append(tv[ind])

    if interface_elements == []:
        interface_elements = np.asarray([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])  # explicit signature required for numba

    interface_elements = np.asarray(interface_elements)

    prn_upd("number of interface elements: {}".format(len(interface_elements)))

    ks_red = np.ones(6 * len(interface_elements))  # reduction of interface shear stiffness for plasticity iteration

    # link interface corner nodes in the same primitive. link0[ind] is connected to link1[ind].
    # This can be done better with a dictionary, but Numba does not support dictionaries.
    link0 = [0]
    link1 = [0]
    for nd in twins:
        for prim in twins[nd]:
            if len(twins[nd][prim]) == 2:
                link0.append(twins[nd][prim][0])
                link1.append(twins[nd][prim][1])

    link0 = np.asarray(link0)
    link1 = np.asarray(link1)

    noce = np.zeros((len(nocoord)), dtype=np.int16)
    # for inno in range(len(nocoord)):
    #     i, j = np.where(elNodes == inno + 1)
    #     noce[inno] = len(i)
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
            bc_ref_nodes = list(set(list(twins.keys())) & set(bcnodes))  # interface reference nodes p/o BC face
            bc_internal_nodes = list(set(bcnodes) - set(list(twins.keys())))  # BC nodes not p/o an interface
            for nd in bc_ref_nodes:
                for prim in twins[nd]:
                    if prim == NiP[bc_internal_nodes[0]]:
                        bcnodes += twins[nd][
                            prim]  # add new interface nodes to the BC face, if they are embedded in the same primitive as the "internal nodes"
            if bcnodes:
                dispfaces.append([bcnodes, bctype, bcvalue])

    dispfaces.append([list(twins.keys()), [False, False, False], [u0, u0, u0]])  # immobilize interface reference nodes

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

    # print(f"fix: {fix}")

    lf = [[0, 0, 0, 0, 0, 0]]  # load face nodes - signature for numba
    pr = [0]  # load face pressure - signature for numba
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
                            if twins:
                                regular_nodes = list(
                                    set(face_nodes) - set(list(twins.keys())))  # nodes outside interface
                                prim = NiP[regular_nodes[0]]  # select primitive of first regular node
                                for index, nd in enumerate(face_nodes):
                                    if nd in twins:
                                        face_nodes[index] = twins[nd][
                                            prim][0]  # swap node number from reference node to new node
                            lf.append(face_nodes)
                            pr.append(sign * obj.Pressure)
                    else:
                        prn_upd("No Faces with Pressure Loads")
    loadfaces = np.asarray(lf)
    pressure = np.array(pr)

    # re-order element nodes
    for el in elNodes:
        el[1], el[2] = el[2], el[1]
        el[4], el[6] = el[6], el[4]
        el[8], el[9] = el[9], el[8]

    return (elNodes, nocoord, fix, fixdof, movdof, loadfaces, materialbyElement, interface_elements, noce, pressure,
            link0, link1, ks_red)


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


# @jit(nopython=True, cache=True)
# def dshp10tet(xi, et, ze, xl):
#     dshp = np.zeros((3, 10), dtype=np.float64)
#     dshpg = np.zeros((3, 10), dtype=np.float64)
#     bmat = np.zeros((6, 30), dtype=np.float64)
#     xs = np.zeros((3, 3), dtype=np.float64)
#     xsi = np.zeros((3, 3), dtype=np.float64)
#
#     # local derivatives of the shape functions: xi-derivative - source: Calculix, G Dhondt
#     dshp[0][0] = 1.0 - 4.0 * (1.0 - xi - et - ze)
#     dshp[0][1] = 4.0 * xi - 1.0
#     # dshp[0][2] = 0.0
#     # dshp[0][3] = 0.0
#     dshp[0][4] = 4.0 * (1.0 - 2.0 * xi - et - ze)
#     dshp[0][5] = 4.0 * et
#     dshp[0][6] = -4.0 * et
#     dshp[0][7] = -4.0 * ze
#     dshp[0][8] = 4.0 * ze
#     # dshp[0][9] = 0.0
#
#     # local derivatives of the shape functions: eta-derivative - source: Calculix, G Dhondt
#     dshp[1][0] = 1.0 - 4.0 * (1.0 - xi - et - ze)
#     # dshp[1][1] = 0.0
#     dshp[1][2] = 4.0 * et - 1.0
#     # dshp[1][3] = 0.0
#     dshp[1][4] = -4.0 * xi
#     dshp[1][5] = 4.0 * xi
#     dshp[1][6] = 4.0 * (1.0 - xi - 2.0 * et - ze)
#     dshp[1][7] = -4.0 * ze
#     # dshp[1][8] = 0.0
#     dshp[1][9] = 4.0 * ze
#
#     # local derivatives of the shape functions: zeta-derivative - source: Calculix, G Dhondt
#     dshp[2][0] = 1.0 - 4.0 * (1.0 - xi - et - ze)
#     # dshp[2][1] = 0.0
#     # dshp[2][2] = 0.0
#     dshp[2][3] = 4.0 * ze - 1.0
#     dshp[2][4] = -4.0 * xi
#     # dshp[2][5] = 0.0
#     dshp[2][6] = -4.0 * et
#     dshp[2][7] = 4.0 * (1.0 - xi - et - 2.0 * ze)
#     dshp[2][8] = 4.0 * xi
#     dshp[2][9] = 4.0 * et
#
#     # xs = np.dot(xl, dshp.T) # local derivative of the global coordinates
#
#     for i in range(3):
#         for j in range(3):
#             xs[i][j] = 0.0
#             for k in range(10):
#                 xs[i][j] += xl[i][k] * dshp[j][k]
#
#     # xsj = np.linalg.det(xs) # Jacobian
#
#     xsj = (xs[0][0] * xs[1][1] * xs[2][2] -
#            xs[0][0] * xs[1][2] * xs[2][1] +
#            xs[0][2] * xs[1][0] * xs[2][1] -
#            xs[0][2] * xs[1][1] * xs[2][0] +
#            xs[0][1] * xs[1][2] * xs[2][0] -
#            xs[0][1] * xs[1][0] * xs[2][2])
#
#     # xsi = np.linalg.inv(xs) # global derivative of the local coordinates
#
#     xsi[0][0] = (xs[1][1] * xs[2][2] - xs[2][1] * xs[1][2]) / xsj
#     xsi[0][1] = (xs[0][2] * xs[2][1] - xs[0][1] * xs[2][2]) / xsj
#     xsi[0][2] = (xs[0][1] * xs[1][2] - xs[0][2] * xs[1][1]) / xsj
#     xsi[1][0] = (xs[1][2] * xs[2][0] - xs[1][0] * xs[2][2]) / xsj
#     xsi[1][1] = (xs[0][0] * xs[2][2] - xs[0][2] * xs[2][0]) / xsj
#     xsi[1][2] = (xs[1][0] * xs[0][2] - xs[0][0] * xs[1][2]) / xsj
#     xsi[2][0] = (xs[1][0] * xs[2][1] - xs[2][0] * xs[1][1]) / xsj
#     xsi[2][1] = (xs[2][0] * xs[0][1] - xs[0][0] * xs[2][1]) / xsj
#     xsi[2][2] = (xs[0][0] * xs[1][1] - xs[1][0] * xs[0][1]) / xsj
#
#     # dshp = np.dot(xsi.T, dshp) # global derivatives of the shape functions
#
#     for i in range(3):
#         for j in range(10):
#             for k in range(3):
#                 dshpg[i][j] += xsi[k][i] * dshp[k][j]
#
#     # computation of the strain interpolation matrix bmat
#     for i in range(10):
#         i3 = 3 * i
#         d00 = dshpg[0][i]
#         d10 = dshpg[1][i]
#         d20 = dshpg[2][i]
#         bmat[0][i3] = d00
#         bmat[1][i3 + 1] = d10
#         bmat[2][i3 + 2] = d20
#         bmat[3][i3] = d10
#         bmat[3][i3 + 1] = d00
#         bmat[4][i3] = d20
#         bmat[4][i3 + 2] = d00
#         bmat[5][i3 + 1] = d20
#         bmat[5][i3 + 2] = d10
#
#     return xsj, dshpg, bmat


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


# @jit(nopython=True, cache=True)
# def hooke(element, materialbyElement):
#     dmat = np.zeros((6, 6), dtype=np.float64)
#     e = materialbyElement[element][0]  # Young's Modulus
#     nu = materialbyElement[element][1]  # Poisson's Ratio
#     dm = e * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu)
#     od = nu / (1.0 - nu)
#     sd = 0.5 * (1.0 - 2.0 * nu) / (1.0 - nu)
#
#     dmat[0][0] = dmat[1][1] = dmat[2][2] = 1.0
#     dmat[3][3] = dmat[4][4] = dmat[5][5] = sd
#     dmat[0][1] = dmat[0][2] = dmat[1][2] = od
#     dmat[1][0] = dmat[2][0] = dmat[2][1] = od
#     dmat *= dm
#
#     return dmat


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


# Nodal point locations
@jit(nopython=True, cache=True)
def nodalPoints():
    # Nodal point locations for a 10-noded tetrahedron
    np10 = np.array([[0.0, 0.0, 0.0],
                     [1.0, 0.0, 0.0],
                     [0.0, 1.0, 0.0],
                     [0.0, 0.0, 1.0],
                     [0.5, 0.0, 0.0],
                     [0.5, 0.5, 0.0],
                     [0.0, 0.5, 0.0],
                     [0.0, 0.0, 0.5],
                     [0.5, 0.0, 0.5],
                     [0.0, 0.5, 0.5]])
    # Nodal point locations for 6-noded triangle + Newton Cotes integration weights
    np6 = np.array([[0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.5, 0.0, 0.166666666666666],
                    [0.5, 0.5, 0.166666666666666],
                    [0.0, 0.5, 0.166666666666666]])
    return np10, np6


# calculate the global stiffness matrix and load vector
@jit(nopython=True, cache=True)
def calcGSM(elNodes, nocoord, materialbyElement, fix, loadfaces, interface_elements,
            grav, kn, ks, pressure, link0, link1, ks_red):
    gp10, gp6 = gaussPoints()
    # np10, np6 = nodalPoints()  # required here only for Newton Coates integration
    ne = len(elNodes)  # number of volume elements
    if interface_elements[0][0] == 0:
        ni = 0
    else:
        ni = len(interface_elements)  # number of interface elements
    nn = len(nocoord[:, 0])  # number of degrees of freedom
    if link0[0] == 0:
        nl = 0
    else:
        nl = len(link0)  # number of link elements

    # number of entries in the lower diagonal of the stiffness matrix: (dof**2 + dof) / 2
    ns = int((30 * 30 + 30) / 2 * ne + (36 * 36 + 36) / 2 * ni + (6 * 6 + 6) / 2 * nl)

    print(f"ne = {ne}")
    print(f"ni = {ni}")
    print(f"nl = {nl}")
    print(f"ns = {ns}\n")

    dmatloc = np.zeros((3, 3), dtype=types.float64)  # local material stiffness matrix
    T = np.zeros((3, 3), dtype=types.float64)  # transformation matrix local-global coordinates
    xlf = np.zeros((3, 6), dtype=types.float64)  # coordinates of load face nodes
    xlv = np.zeros((10, 3), dtype=types.float64)  # coordinates of volume element nodes
    xli = np.zeros((3, 6), dtype=types.float64)  # coordinates of interface element nodes
    nodes_int = np.zeros(12, dtype=types.int64)  # interface node numbers
    # gsm = np.zeros((3 * nn, 3 * nn), dtype=types.float64)  # global stiffness matrix
    glv = np.zeros((3 * nn), dtype=types.float64)  # global load vector
    row = np.zeros(ns, dtype=types.int64)  # row indices of COO matrix
    col = np.zeros(ns, dtype=types.int64)  # column indices of COO matrix
    stm = np.zeros(ns, dtype=types.float64)  # stiffness values of COO matrix
    modf = np.zeros((3 * nn), dtype=types.float64)  # stiffness values of COO matrix
    dof = np.zeros(30, dtype=types.int64)
    dmat = np.zeros((6, 6), dtype=np.float64)
    bmatV = np.zeros((6, 30), dtype=np.float64)
    bmatS = np.zeros((3, 36), dtype=np.float64)

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
            gamma[2::3] += grav * shp * ip[3] * abs(xsj)
            V += xsj * ip[3]  # Element volume - not used

        for i in range(10):
            nd = nodes[i] - 1
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

    # interface stiffness value
    kmax = 0.01 * np.amax(stm)

    # For each interface element calculate the element matrix and
    # add to global stiffness matrix

    for el in range(ni):

        if interface_elements[0][0] == 0: break  # first interface element is a dummy signature for numba
        for i in range(12):
            nodes_int[i] = interface_elements[el][i]
        for i in range(6):
            nd = nodes_int[i]
            for j in range(3):
                xli[j][i] = nocoord[nd - 1][j]

        esm = np.zeros((36, 36), dtype=np.float64)
        dmatloc[0][0] = kn * kmax

        # integrate element matrix (np6: Newton Cotes, gp6: Gauss)

        # ORIGINAL PYTHON
        # for ip in gp6:
        #     xi = ip[0]
        #     et = ip[1]
        #     xsj, shp, bmat, xx, xt, xp = shape6tri(xi, et, xli)
        #     T = np.array([xp, xx, xt])
        #     dmatglob = np.dot(T.T, np.dot(dmatloc, T))
        #     esm += np.dot(bmat.T, np.dot(dmatglob, bmat)) * abs(xsj) * ip[2]

        # NUMBA
        for index in range(len(gp6)):
            dmatloc[1][1] = dmatloc[2][2] = ks_red[6 * el + index] * ks * kmax

            ip = gp6[index]
            xi = ip[0]
            et = ip[1]
            xsj, shp, bmatS, xx, xt, xp = shape6tri(xi, et, xli)
            T[0] = xp
            T[1] = xx
            T[2] = xt
            dmatglob = np.dot(T.T, np.dot(dmatloc,
                                          T))  # transformation of local material stiffness matrix to global coordinates
            esm += np.dot(bmatS.T, np.dot(dmatglob, bmatS)) * abs(xsj) * ip[2]

        # add Element matrix to global stiffness matrix
        # for i in range(12):
        #     iglob = nodes_int[i] - 1
        #     iglob3 = 3 * iglob
        #     i3 = 3 * i
        #     for j in range(12):
        #         jglob = nodes_int[j] - 1
        #         jglob3 = 3 * jglob
        #         j3 = j * 3
        #         for k in range(3):
        #             for l in range(3):
        #                 gsm[iglob3 + k, jglob3 + l] += esm[i3 + k, j3 + l]

        for i in range(12):
            iglob = nodes_int[i] - 1
            i3 = 3 * i
            for j in range(i + 1):
                jglob = nodes_int[j] - 1
                j3 = j * 3
                for k in range(3):
                    for l in range(k + 1):
                        row[pos] = i3
                        col[pos] = j3
                        stm[pos] = esm[i3 + k, j3 + l]
                        pos += 1

    # add stiff links for corner interface nodes in the same primitive
    if (link0[0] != 0):
        stiff = 100 * kmax
        esm = np.array([[stiff, -stiff], [-stiff, stiff]])
        for index, node in enumerate(link0):
            n3a = 3 * (int(node) - 1)
            n3b = 3 * (int(link1[index]) - 1)
            for i in range(2):
                n1 = n3a + i
                for j in range(i + 1):
                    n2 = n3b + j
                    row[pos] = max(n1, n2)
                    col[pos] = min(n1, n2)
                    stm[pos] = esm[i, j]

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

    return stm, row, col, glv, modf, kmax


# calculate load-deflection curve
def calcDisp(elNodes, nocoord, fix, fixdof, movdof, modf, materialbyElement, interface_elements, kmax, stm, row, col,
             glv, nstep, iterat_max, error_max, relax, scale_re, scale_up, scale_dn, sig_yield, phi, psi, ts,
             shr_yield, kn, ks, out_disp, link0, link1, elMat, loadFaces, gravity, pressure, ks_red, mmod):
    ndof = len(glv)  # number of degrees of freedom
    nelem = len(elNodes)  # number of elements
    t0 = time.perf_counter()
    gsm = scsp.csc_matrix((stm, (row, col)), shape=(ndof, ndof))  # construct sparse global stiffness matrix
    t1 = time.perf_counter()
    prn_upd("construct sparse global stiffness matrix: {:.2e} s".format((t1 - t0)))

    # print([list(l) for l in elNodes])
    if interface_elements[0][0] != 0:
        ninter = len(interface_elements)
    else:
        ninter = 0

    # TODO: improve error estimate for pure displacement control
    qnorm = np.linalg.norm(glv)
    if qnorm < 1.0: qnorm = 1.0

    # Cholesky decomposition of the global stiffness matrix and elastic solution
    # using cholmod
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
    sig_old =  np.zeros(24 * nelem, dtype=np.float64)
    peeq_new = np.zeros(4 * nelem, dtype=np.float64)  # equivalent plastic strain in Tet10
    trac_new = np.zeros(18 * max(ninter, 1), dtype=np.float64)  # contact stress in Tri6
    disp_new = np.zeros(ndof, dtype=np.float64)  # displacement results
    lbd = np.zeros(1, dtype=np.float64)  # load level
    rfl = np.zeros(1, dtype=np.float64)  # reaction force level (for displacement control)

    gp10, gp6 = gaussPoints()

    # determine elastic reaction force on moving boundary
    if max(movdof) == 1:
        qelastic = np.zeros(3 * len(nocoord), dtype=np.float64)
        update_stress_load(gp10, elNodes, nocoord, materialbyElement, sig_yield, ue, sig_old, sig_new, peeq_new, qelastic)
        # sig_new, peeq_new, qelastic = \
        #     update_stress_load(gp10, elNodes, nocoord, materialbyElement, sig_yield, ue, sig_new)
        # update_stress_load(elNodes, nocoord, materialbyElement, sig_yield, phi, psi, ts, ue,
        #                     np.array(interface_elements), trac_new, kmax, kn, ks, shr_yield, sig_new,
        #                     peeq_new, link0, link1, ks_red, mmod)

        qelastic *= movdof
        qelnorm = np.linalg.norm(qelastic)
        qnorm = qelnorm
        print("qnorm", qnorm)

    step = -1
    cnt = True

    un = [0.]
    # icall = 0

    if float(nstep) == 1.0:
        # perform an elastic (one-step) analysis
        step = 0
        out_disp = 1
        lbd = np.append(lbd, 1.0)
        rfl = np.append(rfl, 1.0)
        disp_new = ue
        un.append(np.max(np.abs(disp_new)))
        # icall = 1
        cnt = False

    call_time_tot = 0.0
    factor_time_tot = 0.0
    iterat_tot = 0

    while cnt:
        cnt = False
        for _ in (range(nstep)):
            step += 1
            restart = 0
            # disp_old = disp_new.copy()
            sig_old = sig_new.copy()
            peeq_old = peeq_new.copy()
            trac_old = trac_new.copy()

            prn_upd("Step: {}".format(step))
            a = du.copy()  # a: Riks control vector
            lbd = np.append(lbd, lbd[step] + dl)  # lbd: load level

            # update stresses and loads

            #t0 = time.perf_counter()
            qin = np.zeros(3 * len(nocoord), dtype=np.float64)
            update_stress_load(gp10, elNodes, nocoord, materialbyElement, sig_yield, du, sig_old, sig_new, peeq_new, qin)
            # sig_new, peeq_new, qin = \
            #     update_stress_load(gp10, elNodes, nocoord, materialbyElement, sig_yield, du, sig_old)
            # update_stress_load(elNodes, nocoord, materialbyElement, sig_yield, phi, psi, ts, du,
            #                     np.array(interface_elements), trac_old, kmax, kn, ks, shr_yield, sig_old,
            #                     peeq_old, link0, link1, ks_red, mmod)

            #t1 = time.perf_counter()

            # icall += 1
            # call_time_tot += t1 - t0

            # print(f"u_s_l time: {t1-t0}")

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

                # calculate tangent stiffness matrix for Interface Elements
                # if iterat == 1 and interface_elements[0][0] != 0:
                #     prn_upd("Update Tangent Stiffness Matrix for Interface Elements")
                #     stm, row, col, glv, kmax = calcGSM(elNodes, nocoord, elMat, fix, loadFaces,
                #                                        interface_elements,
                #                                        gravity, kn, ks, pressure, link0, link1,
                #                                        ks_red)
                #
                #     gsm = scsp.csc_matrix((stm, (row, col)), shape=(ndof, ndof))
                #
                #     factor = cholesky(gsm)
                #     # K_inv * External Load - required in Riks control
                #     ue = factor(glv)

                # K_inv * Unbalanced Force - displacement correction
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

                # t0 = time.perf_counter()
                qin = np.zeros(3 * len(nocoord), dtype=np.float64)
                update_stress_load(gp10, elNodes, nocoord, materialbyElement, sig_yield, du, sig_old, sig_new, peeq_new, qin)
                # update_stress_load(elNodes, nocoord, materialbyElement, sig_yield, phi, psi, ts, du,
                #                    np.array(interface_elements), trac_old, kmax, kn, ks, shr_yield, sig_old,
                #                    peeq_old, link0, link1, ks_red, mmod)

                # t1 = time.perf_counter()

                # icall += 1
                # call_time_tot += t1 - t0

                # print(f"u_s_l time: {t1 - t0}")

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
                        # du = (disp[step] - disp[step - 1]) / scale_re / restart
                        du = (disp_new - disp_old) / scale_re / restart
                    else:
                        # for first step only
                        dl = dl0 / scale_re / restart
                        du = dl * ue / scale_re / restart
                    lbd[step + 1] = lbd[step] + dl

                    # t0 = time.perf_counter()
                    qin = np.zeros(3 * len(nocoord), dtype=np.float64)
                    update_stress_load(gp10, elNodes, nocoord, materialbyElement, sig_yield, du, sig_old, sig_new, peeq_new, qin)
                    # update_stress_load(elNodes, nocoord, materialbyElement, sig_yield, phi, psi, ts, du,
                    #                    np.array(interface_elements), trac_old, kmax, kn, ks, shr_yield, sig_old,
                    #                    peeq_old, link0, link1, ks_red, mmod)
                    # t1 = time.perf_counter()

                    # icall += 1
                    # call_time_tot += t1 - t0

                    # print(f"u_s_l time: {t1 - t0}")

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

        # plot load-displacement curve - TODO: move to output / post-processing
        if max(movdof) == 1:
            cnt = plot(un, rfl)
        else:
            cnt = plot(un, lbd)

    # prn_upd("call_time_tot: {:.2e}, icall: {}, average time per stress update: {:.2e}".format(call_time_tot, icall,
    #                                                                                           call_time_tot / icall))

    prn_upd("total time evaluating K_inv * r: {}".format(factor_time_tot))
    prn_upd("total number of iterations: {}".format(iterat_tot))
    if iterat_tot != 0:
        prn_upd("average time to evaluate K_inv * r: {}".format(factor_time_tot / iterat_tot))

    # out = min(step + 1, abs(int(out_disp)))
    out = step + 1
    if out_disp > 0:
        u_out = un[out]
        l_out = lbd[out]
        # prn_upd("\n************************************************************\n")
        prn_upd("Step: {0:2d} Load level: {1:.3f} Displacement: {2:.4e}".format(out, l_out,
                                                                                u_out))
        # prn_upd("\n************************************************************\n")
        # return disp[out], sig, peeq, trac
        return disp_new, sig_new, peeq_new, trac_new

    else:
        u_out = un[out] - un[out - 1]
        l_out = lbd[out] - lbd[out - 1]
        # prn_upd("\n************************************************************\n")
        prn_upd("Step: {0:2d} Load level increment: {1:.3f} Displacement increment: {2:.4e}".format(out, l_out,
                                                                                                    u_out))
        # prn_upd("\n************************************************************\n")
        # return disp[out] - disp[out - 1], sig, peeq, trac
        return disp_new - disp_old, sig_new, peeq_new, trac_new


# plot the load-deflection curve
def plot(un, lbd):
    class Index(object):
        def stop(self, event):
            self.cnt = False
            plt.close()
            self.clicked = True

        def add(self, event):
            self.cnt = True
            plt.close()
            self.clicked = True

    callback = Index()
    callback.cnt = False
    callback.clicked = False
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.2)
    ax.plot(un, lbd, '-ok')
    ax.set(xlabel='displacement [mm]', ylabel='load factor [-] or load [N]',
           title='')
    ax.grid()
    axstop = plt.axes([0.7, 0.05, 0.1, 0.075])
    axadd = plt.axes([0.81, 0.05, 0.1, 0.075])
    bstop = Button(axstop, 'stop')
    bstop.on_clicked(callback.stop)
    badd = Button(axadd, 'add')
    badd.on_clicked(callback.add)
    plt.show()

    while True:
        plt.pause(0.01)
        if callback.clicked:
            break

    return callback.cnt


# update stresses and loads

@jit(nopython=True, cache=True, fastmath=True)
# @jit(
#     "types.Tuple((float64[::1], float64[::1], float64[::1]))(float64[:,::1], int64[:,::1], float64[:,::1], float64[:,::1], float64, float64[::1], float64[::1])",
#     nopython=True, cache=True, fastmath=True)
def update_stress_load(gp10, elNodes, nocoord, materialbyElement, sig_yield, du, sig, sig_update, peeq_update, qin):
    u10 = np.empty(30, dtype=np.float64)  # displacements for the 10 tetrahedral nodes
    # sig_update = np.empty(24 * len(elNodes), dtype=np.float64)  # updated stresses
    # qin = np.zeros(3 * len(nocoord), dtype=np.float64)  # internal load vector
    # peeq_update = np.empty(4 * len(elNodes), dtype=np.float64)  # updated peeq
    # eps = np.empty(6, dtype=np.float64)
    sig_test = np.empty(6, dtype=np.float64)
    # bmat = np.zeros((6, 30), dtype=np.float64)
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
    # dshp = np.zeros((3, 10), dtype=np.float64)
    dshp0 = np.empty(10, dtype=np.float64)
    dshp1 = np.empty(10, dtype=np.float64)
    dshp2 = np.empty(10, dtype=np.float64)
    dshpg0 = np.empty(10, dtype=np.float64)
    dshpg1 = np.empty(10, dtype=np.float64)
    dshpg2 = np.empty(10, dtype=np.float64)
    # xs = np.zeros((3, 3), dtype=np.float64)
    # xsi = np.zeros((3, 3), dtype=np.float64)

    # VOLUME ELEMENTS

    hooke(0, materialbyElement, dmat)  # local material stiffness matrix

    for el, nodes in enumerate(elNodes):
        # pass
        # xlv = nocoord[nodes - 1]
        # print(xlv)
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
            # u10[i3:i3+3] = du[n3:n3+3]
        for i in range(4):
            ip = gp10[i]
            ippos = elpos + 6 * i

            # xsj = dshp10tet(ip[0], ip[1], ip[2], xlv, bmat)
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

            # xs = np.dot(xlv, dshp.T) # local derivative of the global coordinates

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

            # for ib in range(3):
            #    for jb in range(3):
            #        tmp = 0.0
            #        for kb in range(10):
            #            tmp += xlv[kb][ib] * dshp[jb][kb]
            #        xs[ib][jb] = tmp

            # xsj = np.linalg.det(xs) # Jacobian

            xsj = (xs00 * xs11 * xs22 -
                   xs00 * xs12 * xs21 +
                   xs02 * xs10 * xs21 -
                   xs02 * xs11 * xs20 +
                   xs01 * xs12 * xs20 -
                   xs01 * xs10 * xs22)

            # xsj = (xs[0][0] * xs[1][1] * xs[2][2] -
            #       xs[0][0] * xs[1][2] * xs[2][1] +
            #       xs[0][2] * xs[1][0] * xs[2][1] -
            #       xs[0][2] * xs[1][1] * xs[2][0] +
            #       xs[0][1] * xs[1][2] * xs[2][0] -
            #       xs[0][1] * xs[1][0] * xs[2][2])

            # xsi = np.linalg.inv(xs) # global derivative of the local coordinates

            xsi00 = (xs11 * xs22 - xs21 * xs12) / xsj
            xsi01 = (xs02 * xs21 - xs01 * xs22) / xsj
            xsi02 = (xs01 * xs12 - xs02 * xs11) / xsj
            xsi10 = (xs12 * xs20 - xs10 * xs22) / xsj
            xsi11 = (xs00 * xs22 - xs02 * xs20) / xsj
            xsi12 = (xs10 * xs02 - xs00 * xs12) / xsj
            xsi20 = (xs10 * xs21 - xs20 * xs11) / xsj
            xsi21 = (xs20 * xs01 - xs00 * xs21) / xsj
            xsi22 = (xs00 * xs11 - xs10 * xs01) / xsj

            # xsi[0][0] = (xs[1][1] * xs[2][2] - xs[2][1] * xs[1][2])
            # xsi[0][1] = (xs[0][2] * xs[2][1] - xs[0][1] * xs[2][2])
            # xsi[0][2] = (xs[0][1] * xs[1][2] - xs[0][2] * xs[1][1])
            # xsi[1][0] = (xs[1][2] * xs[2][0] - xs[1][0] * xs[2][2])
            # xsi[1][1] = (xs[0][0] * xs[2][2] - xs[0][2] * xs[2][0])
            # xsi[1][2] = (xs[1][0] * xs[0][2] - xs[0][0] * xs[1][2])
            # xsi[2][0] = (xs[1][0] * xs[2][1] - xs[2][0] * xs[1][1])
            # xsi[2][1] = (xs[2][0] * xs[0][1] - xs[0][0] * xs[2][1])
            # xsi[2][2] = (xs[0][0] * xs[1][1] - dshpg[1][0] * dshpg[0][1])

            # dshp = np.dot(dshpgi.T, dshp) # global derivatives of the shape functions

            for jb in range(10):
                dshpg0[jb] = xsi00 * dshp0[jb] + xsi10 * dshp1[jb] + xsi20 * dshp2[jb]
                dshpg1[jb] = xsi01 * dshp0[jb] + xsi11 * dshp1[jb] + xsi21 * dshp2[jb]
                dshpg2[jb] = xsi02 * dshp0[jb] + xsi12 * dshp1[jb] + xsi22 * dshp2[jb]

            # for ib in range(3):
            #    for jb in range(10):
            #        tmp = 0.0
            #        for kb in range(3):
            #            tmp += xsi[kb][ib] * dshp[kb][jb]
            #        dshpg[ib][jb] = tmp

            # computation of the strain interpolation matrix bmat
            for ib in range(10):
                i3 = 3 * ib
                d00 = dshpg0[ib]
                d10 = dshpg1[ib]
                d20 = dshpg2[ib]
                # bmat[0][i3] = d00
                # bmat[1][i3 + 1] = d10
                # bmat[2][i3 + 2] = d20
                # bmat[3][i3] = d10
                # bmat[3][i3 + 1] = d00
                # bmat[4][i3] = d20
                # bmat[4][i3 + 2] = d00
                # bmat[5][i3 + 1] = d20
                # bmat[5][i3 + 2] = d10
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

            #            # strain

            eps = np.zeros(6, dtype=np.float64)
            for j in range(10):
                j3 = 3 * j
                eps[0] += bm0[j] * u10[j3]
                eps[1] += bm1[j] * u10[j3 + 1]
                eps[2] += bm2[j] * u10[j3 + 2]
                eps[3] += bm3[j] * u10[j3] + bm4[j] * u10[j3 + 1]
                eps[4] += bm5[j] * u10[j3] + bm6[j] * u10[j3 + 2]
                eps[5] += bm7[j] * u10[j3 + 1] + bm8[j] * u10[j3 + 2]

            # strain
            # for j in range(6):
            #    tmp = 0.0
            #    for k in range(30):
            #        tmp += bmat[j, k] * u10[k]
            #    eps[j] = tmp

            # elastic test stress

            for j in range(6):
                tmp = sig[ippos + j]
                for k in range(6):
                    tmp += dmat[j, k] * eps[k]
                sig_test[j] = tmp

            # sig_update = vmises_original(ippos, sig_test, sig_yield, sig_update)
            sxx, syy, szz, sxy, syz, szx = vmises_original_optimised(sig_test, sig_yield)
            sig_update[ippos:ippos + 6] = sxx, syy, szz, sxy, syz, szx

            ipxsj = ip[3] * abs(xsj)
            for j in range(10):
                j3 = 3 * j
                elv[j3] += (bm0[j] * sxx + bm3[j] * sxy + bm5[j] * syz) * ipxsj
                elv[j3 + 1] += (bm1[j] * syy + bm4[j] * sxy + bm7[j] * szx) * ipxsj
                elv[j3 + 2] += (bm2[j] * szz + bm6[j] * syz + bm8[j] * szx) * ipxsj

            # for j in range(30):
            #    tmp = 0.0
            #    for k in range(6):
            #        tmp += bmat[k, j] * sig_update[ippos + k] * ip[3] * abs(xsj)
            #    elv[j] = tmp
        for i in range(10):
            iglob = nodes[i] - 1
            iglob3 = 3 * iglob
            i3 = 3 * i
            for k in range(3):
                qin[iglob3 + k] += elv[i3 + k]

    return #sig_update, peeq_update, qin


# @jit(nopython=True, cache=True, fastmath=True)
# # @jit(
# #     "types.Tuple((float64[::1], float64[::1], float64[::1]))(float64[:,::1], int64[:,::1], float64[:,::1], float64[:,::1], float64, float64[::1], float64[::1])",
# #     nopython=True, cache=True, fastmath=True)
# def update_stress_load(gp10, elNodes, nocoord, materialbyElement, sig_yield, du, sig):
#     u10 = np.zeros(30, dtype=types.float64)  # displacements for the 10 tetrahedral nodes
#     sig_update = np.zeros(24 * len(elNodes), dtype=types.float64)  # updated stresses
#     qin = np.zeros(3 * len(nocoord), dtype=types.float64)  # internal load vector
#     peeq_update = np.zeros(4 * len(elNodes), dtype=types.float64)  # updated peeq
#     # eps = np.zeros(6, dtype=types.float64)
#     sig_test = np.zeros(6, dtype=types.float64)
#     # bmat = np.zeros((6, 30), dtype=np.float64)
#     bm = np.zeros((9, 10), dtype=np.float64)
#     dmat = np.zeros((6, 6), dtype=np.float64)
#     xlv = np.zeros((10, 3), dtype=np.float64)
#     dshp = np.zeros((3, 10), dtype=np.float64)
#     xs = np.zeros((3, 3), dtype=np.float64)
#     xsi = np.zeros((3, 3), dtype=np.float64)
#     dshpg = np.zeros((3, 10), dtype=np.float64)
#
#     # VOLUME ELEMENTS
#
#     hooke(0, materialbyElement, dmat)  # local material stiffness matrix
#
#     for el, nodes in enumerate(elNodes):
#         # xlv = nocoord[nodes - 1]
#         for i, nd in enumerate(nodes):
#             co = nocoord[nd - 1]
#             xlv[i] = co
#         elpos = 24 * el  # 4 integration points with 6 stress components each
#         elv = np.zeros(30, dtype=np.float64)  # element load vector, 10 nodes with 3 load components each
#
#         for index, nd in enumerate(nodes):
#             n3 = 3 * (nd - 1)
#             i3 = 3 * index
#             u10[i3] = du[n3]
#             u10[i3 + 1] = du[n3 + 1]
#             u10[i3 + 2] = du[n3 + 2]
#         for i in range(4):
#             ip = gp10[i]
#             ippos = elpos + 6 * i
#
#             # xsj = dshp10tet(ip[0], ip[1], ip[2], xlv, bmat)
#
#             # ------------------------------------------------------------------------------------------------------------
#             xi = ip[0]
#             et = ip[1]
#             ze = ip[2]
#             # local derivatives of the shape functions: xi-derivative - source: Calculix, G Dhondt
#             dshp[0][0] = 1.0 - 4.0 * (1.0 - xi - et - ze)
#             dshp[0][1] = 4.0 * xi - 1.0
#             # dshp[0][2] = 0.0
#             # dshp[0][3] = 0.0
#             dshp[0][4] = 4.0 * (1.0 - 2.0 * xi - et - ze)
#             dshp[0][5] = 4.0 * et
#             dshp[0][6] = -4.0 * et
#             dshp[0][7] = -4.0 * ze
#             dshp[0][8] = 4.0 * ze
#             # dshp[0][9] = 0.0
#
#             # local derivatives of the shape functions: eta-derivative - source: Calculix, G Dhondt
#             dshp[1][0] = 1.0 - 4.0 * (1.0 - xi - et - ze)
#             # dshp[1][1] = 0.0
#             dshp[1][2] = 4.0 * et - 1.0
#             # dshp[1][3] = 0.0
#             dshp[1][4] = -4.0 * xi
#             dshp[1][5] = 4.0 * xi
#             dshp[1][6] = 4.0 * (1.0 - xi - 2.0 * et - ze)
#             dshp[1][7] = -4.0 * ze
#             # dshp[1][8] = 0.0
#             dshp[1][9] = 4.0 * ze
#
#             # local derivatives of the shape functions: zeta-derivative - source: Calculix, G Dhondt
#             dshp[2][0] = 1.0 - 4.0 * (1.0 - xi - et - ze)
#             # dshp[2][1] = 0.0
#             # dshp[2][2] = 0.0
#             dshp[2][3] = 4.0 * ze - 1.0
#             dshp[2][4] = -4.0 * xi
#             # dshp[2][5] = 0.0
#             dshp[2][6] = -4.0 * et
#             dshp[2][7] = 4.0 * (1.0 - xi - et - 2.0 * ze)
#             dshp[2][8] = 4.0 * xi
#             dshp[2][9] = 4.0 * et
#
#             # xs = np.dot(xl, dshp.T) # local derivative of the global coordinates
#
#             # for ib in range(3):
#             #     for jb in range(3):
#             #         xs[ib][jb] = 0.0
#             #         for kb in range(10):
#             #             xs[ib][jb] += xlv[kb][ib] * dshp[jb][kb]
#
#             for ib in range(3):
#                 for jb in range(3):
#                     tmp = 0.0
#                     for kb in range(10):
#                         tmp += xlv[kb][ib] * dshp[jb][kb]
#                     xs[ib][jb] = tmp
#
#             # xsj = np.linalg.det(xs) # Jacobian
#
#             xsj = (xs[0][0] * xs[1][1] * xs[2][2] -
#                    xs[0][0] * xs[1][2] * xs[2][1] +
#                    xs[0][2] * xs[1][0] * xs[2][1] -
#                    xs[0][2] * xs[1][1] * xs[2][0] +
#                    xs[0][1] * xs[1][2] * xs[2][0] -
#                    xs[0][1] * xs[1][0] * xs[2][2])
#
#             # xsi = np.linalg.inv(xs) # global derivative of the local coordinates
#
#             xsi[0][0] = (xs[1][1] * xs[2][2] - xs[2][1] * xs[1][2])
#             xsi[0][1] = (xs[0][2] * xs[2][1] - xs[0][1] * xs[2][2])
#             xsi[0][2] = (xs[0][1] * xs[1][2] - xs[0][2] * xs[1][1])
#             xsi[1][0] = (xs[1][2] * xs[2][0] - xs[1][0] * xs[2][2])
#             xsi[1][1] = (xs[0][0] * xs[2][2] - xs[0][2] * xs[2][0])
#             xsi[1][2] = (xs[1][0] * xs[0][2] - xs[0][0] * xs[1][2])
#             xsi[2][0] = (xs[1][0] * xs[2][1] - xs[2][0] * xs[1][1])
#             xsi[2][1] = (xs[2][0] * xs[0][1] - xs[0][0] * xs[2][1])
#             xsi[2][2] = (xs[0][0] * xs[1][1] - xs[1][0] * xs[0][1])
#
#             xsi /= xsj
#
#             # dshp = np.dot(xsi.T, dshp) # global derivatives of the shape functions
#
#             # dshpg = np.zeros((3, 10), dtype=np.float64)
#
#             # for ib in range(3):
#             #     for jb in range(10):
#             #         for kb in range(3):
#             #             dshpg[ib][jb] += xsi[kb][ib] * dshp[kb][jb]
#
#             for ib in range(3):
#                 for jb in range(10):
#                     tmp = 0.0
#                     for kb in range(3):
#                         tmp += xsi[kb][ib] * dshp[kb][jb]
#                     dshpg[ib][jb] = tmp
#
#             # computation of the strain interpolation matrix bmat
#             for ib in range(10):
#                 # i3 = 3 * ib
#                 d00 = dshpg[0][ib]
#                 d10 = dshpg[1][ib]
#                 d20 = dshpg[2][ib]
#                 # bmat[0][i3] = d00
#                 # bmat[1][i3 + 1] = d10
#                 # bmat[2][i3 + 2] = d20
#                 # bmat[3][i3] = d10
#                 # bmat[3][i3 + 1] = d00
#                 # bmat[4][i3] = d20
#                 # bmat[4][i3 + 2] = d00
#                 # bmat[5][i3 + 1] = d20
#                 # bmat[5][i3 + 2] = d10
#                 bm[0][ib] = d00
#                 bm[1][ib] = d10
#                 bm[2][ib] = d20
#                 bm[3][ib] = d10
#                 bm[4][ib] = d00
#                 bm[5][ib] = d20
#                 bm[6][ib] = d00
#                 bm[7][ib] = d20
#                 bm[8][ib] = d10
#
#             # ------------------------------------------------------------------------------------------
#
#             # strain
#             eps = np.zeros(6, dtype=types.float64)
#             for j in range(10):
#                 j3 = 3 * j
#                 eps[0] += bm[0][j] * u10[j3]
#                 eps[1] += bm[1][j] * u10[j3 + 1]
#                 eps[2] += bm[2][j] * u10[j3 + 2]
#                 eps[3] += bm[3][j] * u10[j3] + bm[4][j] * u10[j3 + 1]
#                 eps[4] += bm[5][j] * u10[j3] + bm[6][j] * u10[j3 + 2]
#                 eps[5] += bm[7][j] * u10[j3 + 1] + bm[8][j] * u10[j3 + 2]
#
#             # for j in range(6):
#             #    tmp = 0.0
#             #    for k in range(30):
#             #        tmp += bmat[j, k] * u10[k]
#             #    eps[j] = tmp
#
#             # for j in range(6):
#             #     eps[j] = 0.0
#             #     for k in range(30):
#             #         eps[j] += bmat[j, k] * u10[k]
#
#             # elastic test stress
#             # for j in range(6):
#             #     sig_test[j] = sig[ippos + j]
#             #     for k in range(6):
#             #         sig_test[j] += dmat[j, k] * eps[k]
#
#             for j in range(6):
#                 tmp = sig[ippos + j]
#                 for k in range(6):
#                     tmp += dmat[j, k] * eps[k]
#                 sig_test[j] = tmp
#
#             # sig_update = vmises_original(ippos, sig_test, sig_yield, sig_update)
#
#             sig_update[ippos:ippos + 6] = vmises_original_optimised(sig_test, sig_yield)
#             print(type(vmises_original_optimised(sig_test, sig_yield)))
#
#             ipxsj = ip[3] * abs(xsj)
#             for j in range(10):
#                 j3 = 3 * j
#                 elv[j3] += (bm[0][j] * sig_update[ippos] + bm[3][j] * sig_update[ippos + 3] + bm[5][j] *
#                             sig_update[ippos + 4]) * ipxsj
#                 elv[j3 + 1] += (bm[1][j] * sig_update[ippos + 1] + bm[4][j] * sig_update[ippos + 3] + bm[7][j] *
#                                 sig_update[ippos + 5]) * ipxsj
#                 elv[j3 + 2] += (bm[2][j] * sig_update[ippos + 2] + bm[6][j] * sig_update[ippos + 4] + bm[8][j] *
#                                 sig_update[ippos + 5]) * ipxsj
#
#             # for j in range(30):
#             #     tmp = elv[j]
#             #     for k in range(6):
#             #         tmp += bmat[k, j] * sig_update[ippos + k] * ip[3] * abs(xsj)
#             #     elv[j] = tmp
#
#             # for j in range(30):
#             #     for k in range(6):
#             #         elv[j] += bmat[k, j] * sig_update[ippos + k] * ip[3] * abs(xsj)
#
#         for i in range(10):
#             iglob = nodes[i] - 1
#             iglob3 = 3 * iglob
#             i3 = 3 * i
#             for k in range(3):
#                 qin[iglob3 + k] += elv[i3 + k]
#
#     return sig_update, peeq_update, qin

# @jit(nopython=True, cache=True, fastmath=True)
# def update_stress_load(elNodes, nocoord, materialbyElement, sig_yield, phi, psi, ts, du,
#                        interface_elements, trac, kmax, kn, ks, shr_yield, sig, peeq, link0, link1, ks_red, mmod):
#     gp10, gp6 = gaussPoints()
#
#     # NUMBA:
#     # xlv = np.zeros((3, 10), dtype=types.float64)  # nodal coordinates volume element
#     # nodes = np.zeros(10, dtype=types.int64)  # node numbers single volume element
#     u10 = np.zeros(30, dtype=types.float64)  # displacements for the 10 tetrahedral nodes
#     sig_update = np.zeros(24 * len(elNodes), dtype=types.float64)  # updated stresses
#     trac_update = np.zeros(18 * len(interface_elements), dtype=types.float64)  # updated tractions
#     qin = np.zeros(3 * len(nocoord), dtype=types.float64)  # internal load vector
#     peeq_update = np.zeros(4 * len(elNodes), dtype=types.float64)  # updated peeq
#     eps = np.zeros(6, dtype=types.float64)
#     sig_test = np.zeros(6, dtype=types.float64)
#
#     # VOLUME ELEMENTS
#
#     # for el in prange(len(elNodes)):
#     #     nodes = elNodes[el]
#     for el, nodes in enumerate(elNodes):
#         xlv = nocoord[nodes - 1].T
#         elpos = 24 * el  # 4 integration points with 6 stress components each
#         dmat = hooke(el, materialbyElement)  # local material stiffness matrix
#         elv = np.zeros(30, dtype=np.float64)  # element load vector, 10 nodes with 3 load components each
#
#         for index, nd in enumerate(nodes):
#             n3 = 3 * (nd - 1)
#             i3 = 3 * index
#             u10[i3] = du[n3]
#             u10[i3 + 1] = du[n3 + 1]
#             u10[i3 + 2] = du[n3 + 2]
#         for i in range(4):
#             ip = gp10[i]
#             ippos = elpos + 6 * i
#
#             xsj, dshp, bmat = dshp10tet(ip[0], ip[1], ip[2], xlv)
#
#             # strain
#             for j in range(6):
#                 eps[j] = 0.0
#                 for k in range(30):
#                     eps[j] += bmat[j, k] * u10[k]
#
#             # elastic test stress
#             for j in range(6):
#                 sig_test[j] = sig[ippos + j]
#                 for k in range(6):
#                     sig_test[j] += dmat[j, k] * eps[k]
#
#             # sig_update = vmises_original(ippos, sig_test, sig_yield, sig_update)
#
#             sig_update[ippos:ippos+6] = vmises_original_optimised(sig_test, sig_yield)
#
#             for j in range(30):
#                 for k in range(6):
#                     elv[j] += bmat[k, j] * sig_update[ippos + k] * ip[3] * abs(xsj)
#
#         for i in range(10):
#             iglob = nodes[i] - 1
#             iglob3 = 3 * iglob
#             i3 = 3 * i
#             for k in range(3):
#                 qin[iglob3 + k] += elv[i3 + k]
#
#     return sig_update, peeq_update, trac_update, qin, ks_red


# @jit(nopython=True, cache=True, nogil=True)
# def update_stress_load(elNodes, nocoord, materialbyElement, sig_yield, phi, psi, ts, du,
#                        interface_elements, trac, kmax, kn, ks, shr_yield, sig, peeq, link0, link1, ks_red, mmod):
#     gp10, gp6 = gaussPoints()
#     np10, np6 = nodalPoints()
#
#     # NUMBA:
#     xlv = np.zeros((3, 10), dtype=types.float64)  # nodal coordinates volume element
#     xli = np.zeros((3, 6), dtype=types.float64)  # nodal coordinates interface element
#     dmatloc = np.zeros((3, 3), dtype=types.float64)  # local material stiffness matrix
#     T = np.zeros((3, 3), dtype=types.float64)  # coordinate transformation matrix
#     nodes = np.zeros(10, dtype=types.int64)  # node numbers single volume element
#     nodes_int = np.zeros(12, dtype=types.int64)  # node numbers single interface element
#     u10 = np.zeros((30), dtype=types.float64)  # displacements for the 10 tetrahedral nodes
#     u12 = np.zeros((36), dtype=types.float64)  # displacements for the 12 triangular node pairs
#     sig_update = np.zeros(24 * len(elNodes), dtype=types.float64)  # updated stresses
#     trac_update = np.zeros(18 * len(interface_elements), dtype=types.float64)  # updated tractions
#     qin = np.zeros(3 * len(nocoord), dtype=types.float64)  # internal load vector
#     peeq_update = np.zeros(4 * len(elNodes), dtype=types.float64)  # updated peeq
#
#     # ORIGINAL PYTHON:
#     # xlv = np.zeros((3, 10), dtype=np.float64)
#     # xli = np.zeros((3, 6), dtype=np.float64)
#     # dmatloc = np.zeros((3, 3), dtype=np.float64)
#     # T = np.zeros((3, 3), dtype=np.float64)
#     # nodes = np.zeros(10, dtype=int)
#     # nodes_int = np.zeros(12, dtype=int)
#     # u10 = np.zeros((30), dtype=np.float64)  # displacements for the 10 tetrahedral nodes
#     # u12 = np.zeros((36), dtype=np.float64)  # displacements for the 12 triangular node pairs
#     # sig_update = np.zeros(24 * len(elNodes), dtype=np.float64)
#     # trac_update = np.zeros(18 * len(interface_elements), dtype=np.float64)
#     # qin = np.zeros(3 * len(nocoord), dtype=np.float64)  # internal load vector
#
#     # VOLUME ELEMENTS
#     for el in range(len(elNodes)):
#         for i in range(10):
#             nodes[i] = elNodes[el][i]
#         elpos = 24 * el  # 4 integration points with 6 stress components each
#         eqpos = 4 * el  # 4 integration points with 1 strain component each
#         dmat = hooke(el, materialbyElement)  # local material stiffness matrix
#         G = dmat[3][3]
#         nu = materialbyElement[el][1]
#         bulk = (dmat[0][0] + dmat[0][1] + dmat[0][2]) / 3.0
#         sphi = np.sin(phi)
#         spsi = np.sin(psi)
#         cphi = np.cos(phi)
#         elv = np.zeros((30), dtype=np.float64)  # element load vector, 10 nodes with 3 load components each
#         # xl = np.array([nocoord[nd-1] for nd in nodes]).T
#         for index, nd in enumerate(nodes):
#             for i in range(3):
#                 xlv[i][index] = nocoord[nd - 1][i]
#             n3 = 3 * (nd - 1)
#             i3 = 3 * index
#             u10[i3] = du[n3]
#             u10[i3 + 1] = du[n3 + 1]
#             u10[i3 + 2] = du[n3 + 2]
#         index = -1
#         for i in range(4):
#             index += 1
#             ip = gp10[i]
#             ippos = elpos + 6 * index
#             ippos_eq = eqpos + i
#             xi = ip[0]
#             et = ip[1]
#             ze = ip[2]
#             xsj, dshp, bmat = dshp10tet(xi, et, ze, xlv)
#
#             # elastic test stress
#             sig_test = np.dot(dmat, np.dot(bmat, u10))
#             sig_test[0] += sig[ippos]
#             sig_test[1] += sig[ippos + 1]
#             sig_test[2] += sig[ippos + 2]
#             sig_test[3] += sig[ippos + 3]
#             sig_test[4] += sig[ippos + 4]
#             sig_test[5] += sig[ippos + 5]
#
#             if mmod == 1:  # von Mises / Drucker Prager
#
#                 # sig_update, peeq_update = vmises(ippos, ippos_eq, sig_test, sig_update, peeq, peeq_update, G, bulk,
#                 #                                  sig_yield, sphi, spsi, cphi, ts)
#
#                 sig_update = vmises_original(ippos, sig_test, sig_yield, sig_update)
#
#             if mmod == 2:  # Mohr Coulomb
#                 sig_update = mcoulomb(ippos, sig_test, sig_update, G, bulk, nu,
#                                       sig_yield, sphi, spsi, cphi, ts)
#
#             elv += np.dot(bmat.T, sig_update[ippos:ippos + 6]) * ip[3] * abs(xsj)
#
#         for i in range(10):
#             iglob = nodes[i] - 1
#             iglob3 = 3 * iglob
#             i3 = 3 * i
#             for k in range(3):
#                 qin[iglob3 + k] += elv[i3 + k]
#
#     # INTERFACE ELEMENTS
#
#     for el in range(len(interface_elements)):
#         if interface_elements[0][0] == 0: break
#
#         for i in range(12):
#             nodes_int[i] = interface_elements[el][i]
#
#         elpos = 18 * el  # 6 integration points with 3 stress components each
#
#         # material matrix in local coordinate system
#         dmatloc[0][0] = kn * kmax
#
#         # element nodal displacements
#         for index, nd in enumerate(nodes_int):
#             n3 = 3 * (nd - 1)
#             i3 = 3 * index
#             u12[i3] = du[n3]
#             u12[i3 + 1] = du[n3 + 1]
#             u12[i3 + 2] = du[n3 + 2]
#         inelv = np.zeros((36), dtype=np.float64)  # interface element load vector, 12 nodes with 3 load components each
#
#         for i in range(6):
#             nd = nodes_int[i]
#             for j in range(3):
#                 xli[j][i] = nocoord[nd - 1][j]
#
#         for i in range(6):
#             dmatloc[1][1] = ks_red[6 * el + i] * ks * kmax
#             dmatloc[2][2] = ks_red[6 * el + i] * ks * kmax
#             ip = gp6[i]
#             nppos = elpos + 3 * i
#             xi = ip[0]
#             et = ip[1]
#             xsj, shp, bmat, xx, xt, xp = shape6tri(xi, et, xli)
#
#             T[0] = xp
#             T[1] = xx
#             T[2] = xt
#
#             # initial interface stress in local coordinates
#             int_stress_ini = np.dot(trac[nppos:nppos + 3], T.T)
#
#             # interface strain increment in local coordinates
#             int_strain_inc = np.dot(np.dot(bmat, u12), T.T)
#
#             # test stress increment in local coordinates
#             int_stress_inc = np.dot(dmatloc, int_strain_inc)
#
#             # total test stress in local coordinates
#             int_stress_tot = int_stress_ini + int_stress_inc
#
#             normal = int_stress_tot[0]
#             shear1 = int_stress_tot[1]
#             shear2 = int_stress_tot[2]
#             shear = np.sqrt(shear1 ** 2 + shear2 ** 2)
#             if shear != 0.0:
#                 fac = min(shr_yield / shear, 1.0)
#             else:
#                 fac = 1.0
#
#             if fac < 1.0: ks_red[6 * el + i] = 0.01  # assumed tangent stiffness value
#
#             # updated stress in local coordinates
#             int_stress_loc = np.array([normal, fac * shear1, fac * shear2])
#
#             # updated stress in global coordinates
#             int_stress_glo = np.dot(T.T, int_stress_loc)
#
#             trac_update[nppos:nppos + 3] = int_stress_glo
#
#             # integrate element load vector
#             inelv += np.dot(bmat.T, int_stress_glo) * abs(xsj) * ip[2]
#
#         for i in range(12):
#             iglob = nodes_int[i] - 1
#             iglob3 = 3 * iglob
#             i3 = 3 * i
#             for k in range(3):
#                 qin[iglob3 + k] += inelv[i3 + k]
#
#     # STIFF LINK FORCES
#     stiff = 100 * kmax
#     for index, node in enumerate(link0):
#         if node != 0:
#             n3a = 3 * (int(node) - 1)
#             n3b = 3 * (int(link1[index]) - 1)
#             for i in range(3):
#                 flink = stiff * (du[n3a + i] - du[n3b + i])
#                 qin[n3a + i] += flink
#                 qin[n3b + i] -= flink
#
#     return sig_update, peeq_update, trac_update, qin, ks_red


@jit(nopython=True, cache=True, nogil=True)
def vmises_original(ippos, sig_test, sig_yield, sig_update):
    p = (sig_test[0] + sig_test[1] + sig_test[2]) / 3.0
    sig_mises = np.sqrt(1.5 * ((sig_test[0] - p) ** 2 + (sig_test[1] - p) ** 2 + (sig_test[2] - p) ** 2) +
                        3.0 * (sig_test[3] ** 2 + sig_test[4] ** 2 + sig_test[5] ** 2))

    fac = sig_yield / sig_mises
    if fac > 1.0: fac = 1.0

    sig_update[ippos] = fac * (sig_test[0] - p) + p
    sig_update[ippos + 1] = fac * (sig_test[1] - p) + p
    sig_update[ippos + 2] = fac * (sig_test[2] - p) + p
    sig_update[ippos + 3] = fac * sig_test[3]
    sig_update[ippos + 4] = fac * sig_test[4]
    sig_update[ippos + 5] = fac * sig_test[5]

    return sig_update


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


@jit(nopython=True, cache=True, nogil=True)
def vmises(ippos, ippos_eq, sig_test, sig_update, peeq, peeq_update, G, bulk, sig_yield, sphi, spsi, cphi, ts):
    p = (sig_test[0] + sig_test[1] + sig_test[2]) / 3.0
    sig_mises = np.sqrt(1.5 * ((sig_test[0] - p) ** 2 + (sig_test[1] - p) ** 2 + (sig_test[2] - p) ** 2) +
                        3.0 * (sig_test[3] ** 2 + sig_test[4] ** 2 + sig_test[5] ** 2))

    # DP fit MC at triaxial compression ridge:
    # mf = 6.0 * sphi / (3.0 - sphi)
    # mg = 6.0 * spsi / (3.0 - spsi)
    # cf = 3.0 * sig_yield * cphi / (3.0 - sphi)

    # DP fit MC at triaxial extension ridge:
    mf = 6.0 * sphi / (3.0 + sphi)
    mg = 6.0 * spsi / (3.0 + spsi)
    cf = 3.0 * sig_yield * cphi / (3.0 + sphi)

    # DP averages MC:
    # mf = 18.0 * sphi / (9.0 - sphi ** 2)
    # mg = 18.0 * spsi / (9.0 - spsi ** 2)
    # cf = 9.0 * sig_yield * cphi / (9.0 - sphi ** 2)

    ct = cf - mf * ts
    if ct < 0.0: ct = 0.0

    # radial stress return to yield surface
    fe = sig_mises + mf * p - cf
    if fe < 0.0: fe = 0.0
    dl = fe / (3.0 * G + bulk * mf * mg)
    fac = 1.0 - 3.0 * G * dl / sig_mises
    pressure = (p - dl * bulk * mg)

    if mf * pressure > cf:  # to DP apex
        fac = 0.0
        pressure = cf / mf

    sig_update[ippos] = fac * (sig_test[0] - p) + pressure
    sig_update[ippos + 1] = fac * (sig_test[1] - p) + pressure
    sig_update[ippos + 2] = fac * (sig_test[2] - p) + pressure
    sig_update[ippos + 3] = fac * sig_test[3]
    sig_update[ippos + 4] = fac * sig_test[4]
    sig_update[ippos + 5] = fac * sig_test[5]

    # equivalent plastic strain peeq
    peeq_update[ippos_eq] = peeq[ippos_eq] + max(sig_mises - sig_yield, 0.0) / G / 3.0

    return sig_update, peeq_update


@jit(nopython=True, cache=True, nogil=True)
def mcoulomb(ippos, sig_test, sig_update, G, bulk, nu, sig_yield, sphi, spsi, cphi, ts):
    ps1, ps2, ps3, psdir = calculate_principal_stress(sig_test)

    mu1 = G * (1.0 + sphi * spsi / (1.0 - 2.0 * nu))
    ddgds1 = G * (1.0 + spsi / (1.0 - 2.0 * nu))
    ddgds2 = 2.0 * G * nu * spsi / (1.0 - 2.0 * nu)
    ddgds3 = G * (-1.0 + spsi / (1.0 - 2.0 * nu))

    femc = (ps1 - ps3) + (ps1 + ps3) * sphi - sig_yield * cphi
    femc *= 0.5
    if femc < 0.0: femc = 0.0

    dl = femc / mu1

    ps1r = ps1 - dl * ddgds1
    ps2r = ps2 - dl * ddgds2
    ps3r = ps3 - dl * ddgds3

    # apex = max(ps1r, ps2r, ps3r) * sphi > 0.5 * sig_yield * cphi  # test for apex
    apex = (ps1r + ps3r) * sphi > sig_yield * cphi  # test for apex

    if apex:
        ps1 = ps2 = ps3 = 0.5 * sig_yield * cphi / sphi
    else:
        ridge12 = ps2r > 1.0000001 * ps1r  # test for ps1 = ps2 ridge
        ridge23 = 1.0000001 * ps2r < ps3r  # test for ps2 = ps3 ridge

        if ridge12:
            femctc = (ps2 - ps3) + (ps2 + ps3) * sphi - sig_yield * cphi
            femctc *= 0.5
            mu2 = 0.5 * G * (1.0 - sphi - spsi +
                             sphi * spsi * (1.0 + 2.0 * nu) / (1.0 - 2.0 * nu))
            dl1 = (mu1 * femc - mu2 * femctc) / (
                    mu1 ** 2 - mu2 ** 2)
            dl2 = (mu1 * femctc - mu2 * femc) / (
                    mu1 ** 2 - mu2 ** 2)
            ps1 = ps1 - dl1 * ddgds1 - dl2 * ddgds2
            ps2 = ps2 - dl1 * ddgds2 - dl2 * ddgds1
            ps3 = ps3 - (dl1 + dl2) * ddgds3

        if ridge23:
            femcte = (ps1 - ps2) + (ps1 + ps2) * sphi - sig_yield * cphi
            femcte *= 0.5
            mu2 = 0.5 * G * (1.0 + sphi + spsi +
                             sphi * spsi * (1.0 + 2.0 * nu) / (1.0 - 2.0 * nu))
            dl1 = (mu1 * femc - mu2 * femcte) / (
                    mu1 ** 2 - mu2 ** 2)
            dl2 = (mu1 * femcte - mu2 * femc) / (
                    mu1 ** 2 - mu2 ** 2)
            ps1 = ps1 - (dl1 + dl2) * ddgds1
            ps2 = ps2 - dl1 * ddgds2 - dl2 * ddgds3
            ps3 = ps3 - dl1 * ddgds3 - dl2 * ddgds2

        if not ridge12 and not ridge23:  # regular return
            ps1 = ps1r
            ps2 = ps2r
            ps3 = ps3r

    # f = (ps1 - ps3) + (ps1 + ps3) * sphi - sig_yield * cphi
    # if f > 1.0e-6: print("WARNING: f= ", f)

    sig_update = cartesian(sig_update, ippos, ps1, ps2, ps3, psdir)

    return sig_update


@jit(nopython=True, cache=True, nogil=True)
def calculate_principal_stress(stress_tensor):
    s11 = stress_tensor[0]  # Sxx
    s22 = stress_tensor[1]  # Syy
    s33 = stress_tensor[2]  # Szz
    s12 = stress_tensor[3]  # Sxy
    s31 = stress_tensor[4]  # Szx
    s23 = stress_tensor[5]  # Syz
    sigma = np.array([
        [s11, s12, s31],
        [s12, s22, s23],
        [s31, s23, s33]
    ])

    eigenvalues, eigenvectors = np.linalg.eig(sigma)

    eigenvalues = eigenvalues.real
    eigenvectors = eigenvectors.real

    idx = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    return (eigenvalues[0], eigenvalues[1], eigenvalues[2], eigenvectors)


@jit(nopython=True, cache=True, nogil=True)
def cartesian(sig_update, ippos, ps1, ps2, ps3, psdir):
    sig_update[ippos + 0] = psdir[0, 0] ** 2 * ps1 + psdir[0, 1] ** 2 * ps2 + psdir[0, 2] ** 2 * ps3
    sig_update[ippos + 1] = psdir[1, 0] ** 2 * ps1 + psdir[1, 1] ** 2 * ps2 + psdir[1, 2] ** 2 * ps3
    sig_update[ippos + 2] = psdir[2, 0] ** 2 * ps1 + psdir[2, 1] ** 2 * ps2 + psdir[2, 2] ** 2 * ps3
    sig_update[ippos + 3] = psdir[0, 0] * psdir[1, 0] * ps1 + psdir[0, 1] * psdir[1, 1] * ps2 + psdir[0, 2] * \
                            psdir[1, 2] * ps3
    sig_update[ippos + 4] = psdir[0, 0] * psdir[2, 0] * ps1 + psdir[0, 1] * psdir[2, 1] * ps2 + psdir[0, 2] * \
                            psdir[2, 2] * ps3
    sig_update[ippos + 5] = psdir[1, 0] * psdir[2, 0] * ps1 + psdir[1, 1] * psdir[2, 1] * ps2 + psdir[1, 2] * \
                            psdir[2, 2] * ps3

    return sig_update


# map stresses to nodes
def mapStresses(elNodes, nocoord, interface_elements, disp, sig, peeq, trac, noce):
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
    expm_int = np.zeros((6, 6), dtype=np.float64)  # extrapolation matrix from Integration Points to 6 tri6 nodes
    ipstress = np.zeros((4, 6), dtype=np.float64)  # Tet10 stresses by Gauss point (4 GP and 6 components)
    ippeeq = np.zeros((4, 1), dtype=np.float64)  # Tet10 peeq by Gauss point (4 GP and 1 component)
    iptrac = np.zeros((6, 3), dtype=np.float64)  # Tri6 tractions by integration point (6 GP and 3 components)

    ip10, ip6 = gaussPoints()
    np10, np6 = nodalPoints()

    tet10stress = np.zeros((len(nocoord), 6), dtype=np.float64)
    tet10peeq = np.zeros((len(nocoord)), dtype=np.float64)
    contactpressurevector = np.zeros((len(nocoord), 3), dtype=np.float64)
    contactpressurevalue = np.zeros((len(nocoord)), dtype=np.float64)
    contactshearvector = np.zeros((len(nocoord), 3), dtype=np.float64)

    xp_node = np.zeros((6, 3), dtype=np.float64)  # normal vector in each of the 6 integration points
    xx_node = np.zeros((6, 3), dtype=np.float64)  # Xi tangential vector in each of the 6 integration points
    xt_node = np.zeros((6, 3),
                       dtype=np.float64)  # shear vector perpendicular to the above 2 vectors in each of the 6 integration points

    # step = len(sig) - 1  # last step in the results

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
            # ipstress[index] = sig[step][ippos:ippos + 6]  # ipstress (4x6): 6 stress components for 4 integration points
            # ippeeq[index] = peeq[step][ipposeq]  # ippeeq (4x1): 1 strain component for 4 integration points
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

    # For each interface element map the tractions to element nodes (np6: Newton Cotes, gp6: Gauss)
    # TODO: for extrapolated Gauss point results nodal averaging is required
    for el, nodes in enumerate(interface_elements):
        if interface_elements[0][0] == 0: break
        elpos = 18 * el
        xl = np.array([nocoord[nd - 1] for nd in nodes[:6]]).T
        for index, ip in enumerate(np6):
            xi = ip[0]
            et = ip[1]
            xsj, shp, bmat, xx, xt, xp = shape6tri(xi, et, xl)
            xp_node[index] = xp
            xx_node[index] = xx
            xt_node[index] = xt
            T = np.array([xp, xx, xt])
            ippos = elpos + 3 * index
            # iptrac[index] = np.dot(T, trac[step][ippos:ippos + 3])  # local stresses at integration points
            iptrac[index] = np.dot(T, trac[ippos:ippos + 3])  # local stresses at integration points
            for i in range(6):
                expm_int[index, i] = shp[i]
        expm_int_inv = np.linalg.inv(expm_int)
        nptrac = np.dot(expm_int_inv, iptrac)  # local stresses extrapolated to nodes
        # add local contact stresses to global vectors
        for index, nd in enumerate(nodes[:6]):
            contactpressurevector[nd - 1] = nptrac[index][0] * xp_node[index]
            contactpressurevalue[nd - 1] = nptrac[index][0]
            contactshearvector[nd - 1] = nptrac[index][1] * xx_node[index] + nptrac[index][2] * xt_node[index]

    return tet10stress, tet10peeq, contactpressurevector, contactpressurevalue, contactshearvector


# fill resultobject with results
def pasteResults(doc, elNodes, nocoord, interface_elements, dis, tet10stress, tet10peeq, contactpressurevector,
                 contactpressurevalue, contactshearvector):
    analysis = doc.getObject("Analysis")

    if analysis is None:
        prn_upd("No Analysis object. Please create one first")
        raise SystemExit()

    nn = len(nocoord)  # number of nodes

    resVol = analysis.addObject(ObjectsFem.makeResultMechanical(doc))[0]
    resInt = analysis.addObject(ObjectsFem.makeResultMechanical(doc))[0]

    # VOLUME MESH START
    # volnodes = {}
    # mode_disp_vol = {}
    elements_tetra10 = {}
    mode_results_vol = {}
    results = []

    nodes = range(1, nn + 1)
    coordinates = map(App.Vector, nocoord)
    displacements = map(App.Vector, np.array_split(dis, nn))
    volnodes = dict(zip(nodes, coordinates))
    mode_disp_vol = dict(zip(nodes, displacements))

    # for index, coord in enumerate(nocoord):
    #     # n3 = 3 * index
    #     # volnodes[index + 1] = App.Vector(coord[0], coord[1], coord[2])
    #     volnodes[index + 1] = App.Vector(coord)
    #     # mode_disp_vol[index + 1] = App.Vector(dis[n3], dis[n3 + 1], dis[n3 + 2])
    #     mode_disp_vol[index + 1] = App.Vector(dis[n3:n3 + 3])

    for index, elem in enumerate(elNodes):
        elements_tetra10[index + 1] = (
            elem[0], elem[2], elem[1], elem[3], elem[6], elem[5], elem[4], elem[7], elem[9], elem[8])

    mode_results_vol['disp'] = mode_disp_vol

    # results.append(mode_results_vol)

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

    # numnodes = len(nocoord)
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

    # INTERFACE MESH START
    if interface_elements[0][0] != 0:
        intnodes = {}
        mode_disp_int = {}
        intconnect = {}
        newnode = {}
        oldnode = {}
        elements_tria6 = {}
        mode_results_int = {}
        results = []

        index = 0
        for i, intel in enumerate(interface_elements):
            for nd in intel[:6]:
                if nd not in intconnect:
                    index += 1
                    intconnect[nd] = index
                    intnodes[index] = App.Vector(nocoord[nd - 1][0], nocoord[nd - 1][1], nocoord[nd - 1][2])
                    mode_disp_int[index] = App.Vector(contactpressurevector[nd - 1][0],
                                                      contactpressurevector[nd - 1][1],
                                                      contactpressurevector[nd - 1][2])
                    newnode[nd] = index
                    oldnode[index] = nd

            elements_tria6[i + 1] = (newnode[intel[0]], newnode[intel[1]], newnode[intel[2]], newnode[intel[3]],
                                     newnode[intel[4]], newnode[intel[5]])

        mode_results_int['disp'] = mode_disp_int

        results.append(mode_results_int)

        mint = {
            'Nodes': intnodes,
            'Seg2Elem': {},
            'Seg3Elem': {},
            'Tria3Elem': {},
            'Tria6Elem': elements_tria6,
            'Quad4Elem': {},
            'Quad8Elem': {},
            'Tetra4Elem': {},
            'Tetra10Elem': {},
            'Hexa8Elem': {},
            'Hexa20Elem': {},
            'Penta6Elem': {},
            'Penta15Elem': {},
            'Results': results
        }

        meshint = itf.make_femmesh(mint)

        result_mesh_object_2 = ObjectsFem.makeMeshResult(doc, 'Result_Mesh_Interface')
        result_mesh_object_2.FemMesh = meshint

        resInt.DisplacementVectors = [App.Vector(dis[3 * (oldnode[n + 1] - 1)], dis[3 * (oldnode[n + 1] - 1) + 1],
                                                 dis[3 * (oldnode[n + 1] - 1) + 2]) for n in range(len(oldnode))]
        resInt.DisplacementLengths = [contactpressurevalue[nd - 1] for nd in
                                      oldnode.values()]  # TODO: This is a dirty hack. move contact pressure to its own result object attribute

        resInt.Mesh = result_mesh_object_2
        resInt.NodeNumbers = [int(key) for key in resInt.Mesh.FemMesh.Nodes.keys()]

        resInt = itf.fill_femresult_mechanical(resInt, results)

    # INTERFACE MESH FINISH

    # Add von Mises Stress and Plastic Strain Ratio to the results
    rt.add_von_mises(resVol)
    rt.add_principal_stress_std(resVol)

    doc.recompute()

    return resInt, resVol
