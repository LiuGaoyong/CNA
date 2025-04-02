"""https://v3.pyscal.org/en/latest/examples/14_common_neighbor_analysis.html"""

from ase import Atoms
from pyscal3 import System


def get_cna(atoms: Atoms):
    sys = System(atoms, format="ase")
    sys.analyze.common_neighbor_analysis()
    return sys.atoms.structure
