# -*- coding: utf-8 -*-
import numpy as np
import math as math
import os, sys
from collections import namedtuple
nb_dir = os.path.split(os.getcwd())[0]
if nb_dir not in sys.path:
    sys.path.append(nb_dir)
from energy.parabola import *
from energy.grad_energy import *

# По координатам 4-х точек рассчитывается косинус диэдрального угла, ими образованного
def dihedral(v1, v2, v3, v4):
    v1 = np.array(v1)
    v2 = np.array(v2)
    v3 = np.array(v3)
    v4 = np.array(v4)
    vm = (v2 + v3) / 2
    vm1 = v1 - vm
    vm4 = v4 - vm
    v23 = v3 - v2
    sp1 = vm1 - v23 * (np.dot(vm1, v23) / np.dot(v23, v23))
    sp2 = vm4 - v23 * (np.dot(vm4, v23) / np.dot(v23, v23))
    cos = np.dot(sp1, sp2) / (np.linalg.norm(sp1) * np.linalg.norm(sp2))
    return cos


# По координатам 3-х точек рассчитывается угл, ими образованный
def angle(v1, v2, v3):
    vv1 = np.array(v1) - np.array(v2)
    vv2 = np.array(v3) - np.array(v2)
    cos = np.dot(vv1, vv2) / (np.linalg.norm(vv1) * np.linalg.norm(vv2))
    acos = np.arccos(cos)
    return acos


# Рассчитываются константы f_ij в E_{nonbonded}
def f_ij(neighbours, two_bonds_neigh, three_bonds_neigh, i, j):
    for idx in neighbours[i]:
        if idx == j:
            return 0
    for idx in two_bonds_neigh[i]:
        if idx == j:
            return 0
    for idx in three_bonds_neigh[i]:
        if idx == j:
            return 0.5
    return 1


# Далее вычисляются отдельные составляющие OPLS force field
# Точнее говоря, вычисляются составляющие энергии,
# которые затрагиваются при изменении координат атома с номером num
# В текущей реализации пересчет каждой составляющей энергии
# составляет O(N), где N - число атомов. На самом деле в этом месте
# можно ускороиться по всем слагаемым, кроме E_elst
def E_dihedral(atoms, dihedrals, num, nghb_d):
    # N = len(dihedrals)
    E = 0
    for k in nghb_d:  # range(0,N):
        # if num in dihedrals[k][0:4]:
        cos = dihedral(atoms[dihedrals[k][0]][1:4], atoms[dihedrals[k][1]][1:4],
                       atoms[dihedrals[k][2]][1:4], atoms[dihedrals[k][3]][1:4])
        E = E + dihedrals[k][4] * (1 + cos)
        E = E + 2 * dihedrals[k][5] * (1 - pow(cos, 2))
        E = E + dihedrals[k][6] * (1 + 4 * pow(cos, 3) - 3 * cos)
        E = E + 8 * dihedrals[k][7] * (pow(cos, 2) - pow(cos, 4))
    return E * 4.184


def E_angle(atoms, angles, num, nghb_a):
    # N = len(angles)
    E = 0
    for k in nghb_a:  # range(0,N):
        E = E + angles[k][3] * pow(
            angle(atoms[angles[k][0]][1:4], 
                  atoms[angles[k][1]][1:4], 
                  atoms[angles[k][2]][1:4]
                 ) - angles[k][4],
            2)
    return E * 4.184


def E_bonds(atoms, bonds, num, nghb_b):
    # N = len(bonds)
    E = 0
    for k in nghb_b:  # range(0,N):
        # if num in bonds[k][0:2]:
        r = np.array(atoms[bonds[k][0]][1:4]) - np.array(atoms[bonds[k][1]][1:4])
        E = E + bonds[k][2] * pow(np.linalg.norm(r) - bonds[k][3], 2)
    return E * 4.184


# 1389.38757 – константа, учитывающая 1/(4πε0), множитель 10^(-10)
# – перевод ангстремов в метры,
# заряды электронов и перевод ккал/моль в кДж/моль
def E_elst(atoms, neighbours, two_bonds_neigh, three_bonds_neigh, num):
    N = len(atoms)
    E = 0
    i = num
    for j in range(0, N):
        if j != i:
            f = f_ij(neighbours, two_bonds_neigh, three_bonds_neigh, i, j)
            if f == 0:
                continue
            else:
                r = np.array(atoms[i][1:4]) - np.array(atoms[j][1:4])
                E = E + f * atoms[i][6] * atoms[j][6] / np.linalg.norm(r)
    return E * 1389.38757


# Для получения формулы из Википедии https://en.wikipedia.org/wiki/OPLS:
# Раскрыть скобки и подставить
# A = eps * (sig ^ 6)
# C = eps * (sig ^ 3)
# Формула переделана для уменьшения числа арифметических операций
# и соответствия параметрам, которые предоставляет OPLS
# 4 * 4.184 – перевод ккал/моль в кДж/моль
def E_vdw(atoms, neighbours, two_bonds_neigh, three_bonds_neigh, num):
    N = len(atoms)
    E = 0
    i = num
    #for j in range(0, N):
    for j in range(i + 1, N):
        if j != i:
            f = f_ij(neighbours, two_bonds_neigh, three_bonds_neigh, i, j)
            if f == 0:
                continue
            else:
                eps = math.sqrt(atoms[i][5] * atoms[j][5])
                sig = atoms[i][4] * atoms[j][4]
                dist = np.linalg.norm(np.array(atoms[i][1:4]) - np.array(atoms[j][1:4]))
                sigSqDivR = (sig / dist ** 2) ** 3
                E = E + f * eps * sigSqDivR * (sigSqDivR - 1)
    return E * 4 * 4.184


def E(atoms, dihedrals, angles, bonds, neighbours, two_bonds_neigh, three_bonds_neigh, num, nghb_d, nghb_a, nghb_b, nghb_nb):
    eng = namedtuple("eng", ["E", "dih", "ang", "bond", "elst", "vdw"])
    E_dih = E_dihedral(atoms, dihedrals, num, nghb_d)
    E_ang = E_angle(atoms, angles, num, nghb_a)
    E_bond = E_bonds(atoms, bonds, num, nghb_b)
    E_elst_ = E_elst(atoms, neighbours, two_bonds_neigh, three_bonds_neigh, num)
    E_vdw_ = E_vdw(atoms, neighbours, two_bonds_neigh, three_bonds_neigh, num)
    E = E_dih + E_ang + E_bond + E_elst_ + E_vdw_
    return eng(E, E_dih, E_ang, E_bond, E_elst_, E_vdw_)


def nghb_dihedrals(dihedrals, num):
    nghb_d = []
    for k in range(0, len(dihedrals)):
        if num in dihedrals[k][0:4]:
            nghb_d.append(k)
    return nghb_d


def nghb_angles(angles, num):
    nghb_a = []
    for k in range(0, len(angles)):
        if num in angles[k][0:3]:
            nghb_a.append(k)
    return nghb_a


def nghb_bonds(bonds, num):
    nghb_b = []
    for k in range(0, len(bonds)):
        if num in bonds[k][0:2]:
            nghb_b.append(k)
    return nghb_b


def nghbnb(atoms, num, r):
    nghb = []
    nghbc = []
    # Список nghb формируется из соседей атома num, находящихся на расстояние менее r
    for k in range(0, len(atoms)):
        if np.linalg.norm(np.array(atoms[k][1:4]) - np.array(atoms[num][1:4])) < r:
            nghb.append(k)
        else:
            nghbc.append(k)
    return nghb, nghbc


def coef_nb(atoms, num, nghbc):
    E = np.array([0.0, 0.0, 0.0])
    i = num
    for j in nghbc:
        r = np.array(atoms[j][1:4]) - np.array(atoms[i][1:4])
        E[0] = E[0] + atoms[i][6] * atoms[j][6] * (atoms[j][1] - atoms[i][1]) / pow(np.linalg.norm(r), 3)
        E[1] = E[1] + atoms[i][6] * atoms[j][6] * (atoms[j][2] - atoms[i][2]) / pow(np.linalg.norm(r), 3)
        E[2] = E[2] + atoms[i][6] * atoms[j][6] * (atoms[j][3] - atoms[i][3]) / pow(np.linalg.norm(r), 3)
    for enj in range(len(E)):
        if math.isnan(E[enj]):
            E[enj] = 0
    return E * 1389.38757