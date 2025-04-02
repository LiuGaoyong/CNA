"""https://asap3.readthedocs.io/en/latest/examples/Running_Common_Neighbor_Analysis.html

Result follwed this paper: https://doi.org/10.1016/0927-0256(94)90109-0
"""

from asap3.analysis.localstructure import FullCNA
from ase import Atoms


def get_cna(atoms: Atoms):
    cna = FullCNA(atoms)
    return cna.get_normal_cna()
