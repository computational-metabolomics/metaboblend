import copy
import itertools
import networkx as nx
import numpy
from operator import itemgetter
from rdkit import Chem


def subset_sum(l, mass, toll=0.001):
    if mass < -toll:
        return
    elif len(l) == 0:
        if -toll <= mass <= toll:
            yield []
        return

    elif abs(sum(l) - mass) <= toll:
        # print(abs(sum(l) - mass), sum(l) - mass)
        yield l
        return

    for subset in subset_sum(l[1:], mass):
        yield subset
    for subset in subset_sum(l[1:], mass - l[0]):
        yield [l[0]] + subset


def combine_ecs(ss2_grp, heavy_atoms, db, accuracy=None, ppm=None):
    ecs = []

    for i in range(len(ss2_grp)):
        if ppm is None:
            atoms = db.select_ecs(ss2_grp[i], heavy_atoms, accuracy)
        else:
            atoms = db.select_ecs(ss2_grp[i], heavy_atoms, accuracy, ppm=ppm)

        if len(atoms) == 0:
            return []

        ecs.append(atoms)

    return ecs


def reindex_atoms(records):
    atoms_available, atoms_to_remove, bond_types = [], [], {}

    mol_comb = Chem.Mol()
    index_atoms = []
    c = 0
    for record in records:
        idxs = []
        for atom in record["mol"].GetAtoms():

            newIdx = atom.GetIdx() + c
            idxs.append(newIdx)

            if atom.GetIdx() in record["degree_atoms"]:
                atoms_available.append(newIdx)
            if atom.GetIdx() in record["dummies"]:
                atoms_to_remove.append(newIdx)
            if atom.GetIdx() in record["bond_types"]:
                bond_types[newIdx] = record["bond_types"][atom.GetIdx()]

        mol_comb = Chem.rdmolops.CombineMols(mol_comb, record["mol"])
        index_atoms.append(idxs)
        c = idxs[-1] + 1

    return mol_comb, atoms_available, atoms_to_remove, bond_types


def add_bonds(mols, edges, atoms_available, bond_types, debug=False):
    rdkit_bond_types = {1: Chem.rdchem.BondType.SINGLE,
                        1.5: Chem.rdchem.BondType.AROMATIC,
                        2: Chem.rdchem.BondType.DOUBLE}

    G = nx.Graph()
    G.add_edges_from(edges)

    if debug:
        print("## Edges from isomorphism:", edges)
        print("## Matching:", sorted(G.nodes()), atoms_available)

    G = nx.relabel_nodes(G, dict(zip(sorted(G.nodes()), atoms_available)))

    # Draw.MolToFile(mols, "main_before_" + "-".join(map(str, atoms_available)) + '.png')
    # Draw.MolToFile(mols, "main_before_indexed_" + "-".join(map(str, atoms_available)) + '.png', includeAtomNumbers=True)

    mol_edit = Chem.EditableMol(mols)
    for edge in G.edges():

        if edge[0] in bond_types:
            bt_start = copy.copy(bond_types[edge[0]])
        else:
            if debug:
                print("## Nested dummy with index: {} ".format(edge[0]))
                print("")
            return None

        if edge[1] in bond_types:
            bt_end = copy.copy(bond_types[edge[1]])
        else:
            if debug:
                print("## Nested dummy with index: {} ".format(edge[1]))
                print("")
            return None

        bondMatches = list(set(bt_start).intersection(bt_end))

        if len(bondMatches) == 0:
            if debug:
                print("## bondMatches empty")
                print("")
            return None
            # yield None
        else:
            bt_start.remove(bondMatches[0])
            bt_end.remove(bondMatches[0])

        # print(edge[0], edge[1], rdkit_bond_types[bondMatches[0]])
        try:
            mol_edit.AddBond(edge[0], edge[1], rdkit_bond_types[bondMatches[0]])
        except KeyError:
            if debug:
                print("Unknown bond type")
            return None

    return mol_edit


def standard_build(mc, exact_mass, db, fn_out, heavy_atoms, max_valence, accuracy, debug=False):

    exact_mass__1 = round(exact_mass)
    exact_mass__0_0001 = round(exact_mass, 4)

    mass_values = db.select_mass_values(str(accuracy), heavy_atoms, max_valence, [])
    subsets = list(subset_sum(mass_values, exact_mass__1))

    configs_iso = db.k_configs()
    out = open(fn_out, "w")

    if debug:
        print("First round (mass: {}) - Values: {} - Correct Sums: {}".format(exact_mass__1, len(mass_values), len(subsets)))
        print("------------------------------------------------------")
    for ss_grp in subsets:

        mass_values_r2 = db.select_mass_values("0_0001", heavy_atoms, max_valence, ss_grp)
        subsets_r2 = list(subset_sum(mass_values_r2, exact_mass__0_0001))

        if debug:
            print("Second round (mass: {}) - Values: {} - Correct Sums: {}".format(exact_mass__0_0001, len(mass_values_r2), len(subsets_r2)))
            print("------------------------------------------------------")

        build_from_subsets(configs_iso, subsets_r2, mc, db, out, heavy_atoms, debug)

    out.close()
    
    
def prescribed_build(mc, exact_mass, db, fn_out, heavy_atoms, max_valence, accuracy, fragment_mass, ppm, debug=False):
    loss = exact_mass - fragment_mass
    exact_mass__1 = round(loss)
    exact_mass__0_0001 = round(loss, 4)

    tolerance = (loss / 1000000) * ppm
    if tolerance < 0.001:
        tolerance = 0.001
    else:
        tolerance = round(tolerance, 4)

    mass_values = db.select_mass_values(str(accuracy), heavy_atoms, max_valence, [])
    subsets = list(subset_sum(mass_values, exact_mass__1))

    configs_iso = db.k_configs()
    out = open(fn_out, "w")

    for ss_grp in subsets:

        mass_values_r2 = db.select_mass_values("0_0001", heavy_atoms, max_valence, ss_grp)
        subsets_r2 = list(subset_sum(mass_values_r2, exact_mass__0_0001, tolerance))

        # append fragment to subsets
        for i, subset in enumerate(subsets_r2):
            subsets_r2[i] = [round(exact_mass - loss, 4)] + subset

        build_from_subsets(configs_iso, subsets_r2, mc, db, out, heavy_atoms, ppm, debug)


def build_from_subsets(configs_iso, subsets_r2, mc, db, out, heavy_atoms, ppm=None, debug=False):
    for ss2_grp in subsets_r2:
        if ppm is None:
            list_ecs = combine_ecs(ss2_grp, heavy_atoms, db, "0_0001")
        else:
            list_ecs = combine_ecs(ss2_grp, heavy_atoms, db, "0_0001", ppm)

        if len(list_ecs) == 0:
            continue

        iii = 0
        for l in itertools.product(*list_ecs):

            sum_ec = list(numpy.array(l).sum(axis=0))
            iii += 1

            if mc != sum_ec and debug:
                print("No match for elemental composition: {}".format(str(sum_ec)))

            elif mc == sum_ec:

                if debug:
                    print("Match elemental composition: {}".format(str(sum_ec)))

                ll = db.select_sub_structures(l)

                if len(ll) == 0:
                    if debug:
                        print("## No substructures found")
                    continue
                elif len(ll) == 1:
                    if debug:
                        print("## Single substructure")
                else:
                    if debug:
                        print("## {} {} substructures found".format(sum([len(subs) for subs in ll]),
                                                                    str([len(subs) for subs in ll])))
                    # print([[sub[""smiles] for sub in subs] for subs in ll])

                if debug:
                    print("## {} substructure combinations".format(len(list(itertools.product(*ll)))))

                for lll in itertools.product(*ll):

                    if debug:
                        for record in lll:
                            print(record)
                        print("---------------")

                    lll = sorted(lll, key=itemgetter('atoms_available', 'valence'))
                    nA, v, vA = (), (), ()
                    for d in lll:
                        nA = nA + (d["atoms_available"],)
                        v = v + (d["valence"],)
                        vA = vA + (tuple(d["degree_atoms"].values()),)

                    if debug:
                        print(str(vA))
                        print("============")
                    # print(configs_iso)
                    # print("============")

                    if str(vA) not in configs_iso:
                        if debug:
                            print("NO:", (str(nA), str(v), str(vA)))
                        continue
                    else:
                        if debug:
                            print("YES:", (str(nA), str(v), str(vA)))

                    # print("## ConnectivityGraphs found (%s)" % (len(list(db.isomorphismGraphs(str(tuple(nA)), str(tuple(v)))))))
                    # print("## Atoms available (n) %s / Valence %s" % (str(tuple(nA)), str(tuple(v))))

                    mol_comb, atoms_available, atoms_to_remove, bond_types = reindex_atoms(lll)
                    if debug:
                        print("## Mols (in memory):", mol_comb)
                        print("## Atoms Available (indexes):", atoms_available)
                        print("## Atoms to remove (dummies):", atoms_to_remove)
                        print("## Type of bonds to form:", bond_types)
                    iso_n = 0
                    for edges in db.isomorphism_graphs(configs_iso[str(vA)]):  # EDGES

                        iso_n += 1
                        if debug:
                            print("## ISO {}".format(iso_n))

                        if debug:
                            print(edges)
                            print("1: Add bonds")
                        mol_e = add_bonds(mol_comb, edges, atoms_available, bond_types)
                        if mol_e is None:
                            continue
                        if debug:
                            print("2: Add bonds")

                        atoms_to_remove.sort(reverse=True)
                        [mol_e.RemoveAtom(a) for a in atoms_to_remove]

                        molOut = mol_e.GetMol()
                        try:
                            Chem.Kekulize(molOut)
                        except:
                            if debug:
                                print("Can't kekulize mol ISO: {}".format(iso_n))
                            continue

                        # Draw.MolToFile(molOut, "main_after_" + "-".join(map(str, atoms_available)) + '.png')
                        try:
                            out.write("{}\t{}\n".format(Chem.MolToSmiles(molOut, kekuleSmiles=True),
                                                        str([item["smiles"] for item in lll])))
                        except RuntimeError:
                            if debug:
                                print("Bad bond type violation")
                        if debug:
                            print("## smi (result): {}".format(
                                Chem.MolToSmiles(molOut, kekuleSmiles=True)))  # , bond_types
