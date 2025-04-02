"""https://github.com/tuy15/Common-Neighbor-Analysis/tree/master"""

from pathlib import Path

import numpy as np
from ase import Atoms


def calculateDistance(a, b):
    return np.linalg.norm(a - b)


def CommonNeighbor(a, b):
    if a[0] != b[0] and calculateDistance(a, b) < 3.5:  # 3.3 for cubo
        return int(b[0] - 1)
    else:
        return 1000000


def Generatelist(a, List):
    NN = []
    for b in List:
        if CommonNeighbor(a, b) != 1000000:
            NN.append(CommonNeighbor(a, b))
    return NN


def Listneighbor(List):
    LN = []
    for a in range(0, len(List)):
        for b in range(0, len(List)):
            if (
                a != b and CommonNeighbor(List[a], List[b]) != 1000000
            ):  # and CommonNeighbor(a,b)!=1000000):
                LN.append([a, b])
    return LN


def get_cna(atoms: Atoms):
    fname = Path(__file__).parent.joinpath("tuy15.txt")
    INDEX = np.genfromtxt(fname.as_posix(), dtype=int).tolist()
    data = atoms.positions
    TYPE = []
    for x in range(0, len(data)):
        Allneighbor = Generatelist(data[x], data)
        # print('NN')  #print(x,len(Allneighbor),Allneighbor)
        Allneighbor_list = [data[i] for i in Allneighbor]
        # print(Allneighbor_list)   # print('CommonNN')
        L = [0, 0, 1, 2, 1, 5, 5]
        LEN = []
        LEN1 = []
        LEN2 = []
        for y in range(0, len(Allneighbor_list)):
            Nextneighbor = Generatelist(Allneighbor_list[y], Allneighbor_list)
            LEN.append(len(Nextneighbor))
            Nextneighbor_list = [data[i] for i in Nextneighbor]
            bonds = Listneighbor(Nextneighbor_list)
            preL = len(np.unique(bonds))  # print(Nextneighbor_list)
            LEN1.append(len(bonds) // 2)
            LEN2.append(L[preL])  # print(bonds,np.unique(bonds));
        # g.write(str(x)+'  '+str(len(Allneighbor))+'  '+str(max(LEN2))+'  '+str(min(LEN2))+'\n')
        G = [len(Allneighbor), max(LEN2), min(LEN2)]
        if G in INDEX:
            TYPE.append(INDEX.index(G))
        else:
            TYPE.append(len(INDEX))
        print(x, G)

    return TYPE
