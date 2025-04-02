"""https://freud.readthedocs.io/en/stable/gettingstarted/examples/examples/NetworkX-CNA.html


Result followed this paper: https://pubs.acs.org/doi/10.1021/j100303a014
"""

from collections import Counter, defaultdict

import networkx as nx
from asap3.analysis.localstructure import GuessLatticeConstant
from ase import Atoms
from freud import AABBQuery, Box


def get_cna(atoms: Atoms):
    box = Box.from_matrix(atoms.cell.array)
    points = atoms.positions
    aq = AABBQuery(box, points)
    rmax = GuessLatticeConstant(atoms) * 0.825
    nl = aq.query(
        points,
        {
            "num_neighbors": 12,
            "exclude_ii": True,
            "r_max": rmax,
        },
    ).toNeighborList()

    # Get all sets of common neighbors.
    common_neighbors = defaultdict(list)
    for i, p in enumerate(points):
        for j in nl.point_indices[nl.query_point_indices == i]:
            for k in nl.point_indices[nl.query_point_indices == j]:
                if i != k:
                    common_neighbors[(i, k)].append(j)

    diagrams = defaultdict(list)
    particle_counts = defaultdict(Counter)
    for (a, b), neighbors in common_neighbors.items():
        # Build up the graph of connections between the
        # common neighbors of a and b.
        g = nx.Graph()
        for i in neighbors:
            for j in set(nl.point_indices[nl.query_point_indices == i]).intersection(
                neighbors
            ):
                g.add_edge(i, j)
        # Define the identifiers for a CNA diagram:
        # The first integer is 1 if the particles are bonded, otherwise 2
        # The second integer is the number of shared neighbors
        # The third integer is the number of bonds among shared neighbors``
        # The fourth integer is an index, just to ensure uniqueness of diagrams
        diagram_type = 2 - int(b in nl.point_indices[nl.query_point_indices == a])
        key = (diagram_type, len(neighbors), g.number_of_edges())
        # If we've seen any neighborhood graphs with this signature,
        # we explicitly check if the two graphs are identical to
        # determine whether to save this one. Otherwise, we add
        # the new graph immediately.
        if key in diagrams:
            isomorphs = [nx.is_isomorphic(g, h) for h in diagrams[key]]
            if any(isomorphs):
                idx = isomorphs.index(True)
            else:
                diagrams[key].append(g)
                idx = diagrams[key].index(g)
        else:
            diagrams[key].append(g)
            idx = diagrams[key].index(g)
        cna_signature = key + (idx,)
        particle_counts[a].update([cna_signature])

    return particle_counts
