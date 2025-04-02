from pprint import pprint

import pytest
from ase import Atoms
from ase.cluster import Octahedron

from cna.asap import get_cna as get_cna_by_asap3
from cna.freud import get_cna as get_cna_by_freud
from cna.pyscal3 import get_cna as get_cna_by_pyscal3


@pytest.mark.parametrize("f", [get_cna_by_asap3, get_cna_by_freud, get_cna_by_pyscal3])
def test_oct(f):
    atoms = Octahedron("Cu", 3)
    atoms = Atoms(atoms.numbers, atoms.positions, cell=[100, 100, 100], pbc=True)
    print(len(atoms))
    pprint(f(atoms))
    print()
