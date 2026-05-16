"""SMILES → NetworkX graph (hydrogen-suppressed, simple unweighted).

Convention: nodes = heavy atoms, edges = bonds (bond order ignored).
This is the standard 'molecular graph' used for classical topological indices.
"""
from __future__ import annotations
import networkx as nx
from rdkit import Chem


def smiles_to_graph(smiles: str) -> nx.Graph | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # Drop counterions / multi-component SMILES: keep largest fragment.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())

    G = nx.Graph()
    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(), symbol=atom.GetSymbol())
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    G = nx.convert_node_labels_to_integers(G)

    if G.number_of_nodes() < 2 or G.number_of_edges() == 0:
        return None
    if not nx.is_connected(G):
        return None
    return G


def _test():
    cases = {
        "C": None,                       # methane: single heavy atom -> None
        "CC": (2, 1),                    # ethane
        "CCO": (3, 2),                   # ethanol
        "c1ccccc1": (6, 6),              # benzene
        "CC(C)C": (4, 3),                # isobutane
        "[Na+].[Cl-]": None,             # salt: largest frag is 1 atom -> None
        "INVALID": None,
    }
    for smi, expected in cases.items():
        G = smiles_to_graph(smi)
        if G is None:
            got = None
        else:
            got = (G.number_of_nodes(), G.number_of_edges())
        assert got == expected, f"{smi}: got {got}, expected {expected}"
        print(f"  {smi!r:18s} -> {got}")
    print("mol_to_graph self-test passed.")


if __name__ == "__main__":
    _test()
