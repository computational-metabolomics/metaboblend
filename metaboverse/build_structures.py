#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019-2020 Ralf Weber
#
# This file is part of MetaboVerse.
#
# MetaboVerse is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MetaboVerse is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MetaboVerse.  If not, see <https://www.gnu.org/licenses/>.
#

import os
import multiprocessing
import copy
import itertools
from functools import partial
import networkx as nx
import numpy
from operator import itemgetter
from rdkit import Chem
from .databases import SubstructureDb


def subset_sum(l, mass, toll=0.001):
    if mass < -toll:
        return

    elif len(l) == 0:
        if -toll <= mass <= toll:
            yield []
        return

    elif abs(sum(l) - mass) <= toll:
        yield l
        return

    for subset in subset_sum(l[1:], mass):
        yield subset
    for subset in subset_sum(l[1:], mass - l[0]):
        yield [l[0]] + subset


def combine_ecs(ss2_grp, db, table_name, accuracy, ppm=None):
    ecs = []

    for i in range(len(ss2_grp)):
        atoms = db.select_ecs(ss2_grp[i], table_name, accuracy, ppm=ppm)

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

        else:
            bt_start.remove(bondMatches[0])
            bt_end.remove(bondMatches[0])

        try:
            mol_edit.AddBond(edge[0], edge[1], rdkit_bond_types[bondMatches[0]])
        except KeyError:
            if debug:
                print("Unknown bond type")
            return None

    return mol_edit


def build(mc, exact_mass, fn_out, heavy_atoms, max_valence, accuracy, max_atoms_available, max_n_substructures,
          path_db_k_graphs="../databases/k_graphs.sqlite", path_pkls="../databases/pkls",
          path_db="../databases/substructures.sqlite", fragment_mass=None, ppm=None, debug=False, out_mode="w"):
    """
    Workflow for generating molecules of a given mass using substructures and connectivity graphs. Can optionally
    take a "prescribed" fragment mass to further filter results; this can be used to incorporate MSn data. Final
    molecules are written to the specified file in smiles format, which the substructures utilised to generate them.

    :param mc: List of integers detailing the molecular composition of the target metabolite, in the format
        [C, H, N, O, P, S].

    :param exact_mass: The exact mass of the target metabolite.

    :param fn_out: The path of the file to which smile strings should be written representing the final structures
        generated; the substructures used to generate these final structures are also written.

    :param heavy_atoms: A list of integers containing the heavy atom counts to consider in substructures to be used
        to build final structures.

    :param max_valence: The maximal total bond orders of substructures to be considered to build final structures
        (ie, the product of `atoms_available` and the degree of their bonds).

    :param max_atoms_available: The maximal atoms available of substructures to be considered for building molecules.

    :param accuracy: To which decimal places of accuracy results are to be limited to. This is only relevant for the
        intial subset_sum pass; in the second stage, the maximum accuracy (d.p.) is always used.

            * **1** Integer level

            * **0_1** One decimal place

            * **0_01** Two decimal places

            * **0_001** Three decimal places

            * **0_0001** Four decimal places

    :param max_n_substructures: The maximum number of substructures to be used for building molecules.

    :param path_db: The path to the SQLite 3 substructure database, as generated by
        :py:meth:`metaboverse.databases.SubstructureDb`.

    :param path_db_k_graphs: The path to the SQLite 3 connectivity database, as generated by
        :py:meth:`metaboverse.databases.create_isomorphism_database`.

    :param path_pkls: The path to the connectivity graphs described by the SQLite 3 connectivity database, as generated
        by :py:meth:`metaboverse.databases.create_isomorphism_database`.

    :param fragment_mass: A mass by which to filter results; if not provided, all possible structures will be generated.

    :param ppm: The allowable error of the query (in parts per million). Designed to be used in accordance with
        `fragment_mass`.

    :param out_mode: The mode in which to write to the output file.

        * **"w"** Create a new file, overwriting any existing file.

        * **"w"** Append results to an existing file.

    :param debug: Debug print statements provide further information on how the function is generating the connectivity
        database.
        * **True** Print debug statements.
        * **False** Hide debug print statements.
    """

    db = SubstructureDb(path_db, path_pkls, path_db_k_graphs)
    table_name = gen_subs_table(db, heavy_atoms, max_valence, max_atoms_available)

    if fragment_mass is None:  # standard build method
        exact_mass__1 = round(exact_mass)
        exact_mass__0_0001 = round(exact_mass, 4)

        tolerance = 0.001
    else:  # prescribed substructure build method
        loss = exact_mass - fragment_mass
        exact_mass__1 = round(loss)
        exact_mass__0_0001 = round(loss, 4)

        tolerance = (loss / 1000000) * ppm
        if tolerance < 0.001:
            tolerance = 0.001
        else:
            tolerance = round(tolerance, 4)

        max_n_substructures -= 1

    if os.name == "nt":  # multiprocessing freeze support on windows
        multiprocessing.freeze_support()

    if out_mode == "w":
        open(fn_out, "w").close()

    mass_values = db.select_mass_values(str(accuracy), [], table_name)

    subsets = list(subset_sum(mass_values, exact_mass__1))

    configs_iso = db.k_configs()
    out = open(fn_out, out_mode)

    if debug:
        print("First round (mass: {}) - Values: {} - Correct Sums: {}".format(exact_mass__1, len(mass_values),
                                                                              len(subsets)))
        print("------------------------------------------------------")
    for ss_grp in subsets:
        if len(ss_grp) > max_n_substructures:
            continue

        mass_values_r2 = db.select_mass_values("0_0001", ss_grp, table_name)
        subsets_r2 = list(subset_sum(mass_values_r2, exact_mass__0_0001, tolerance))

        if fragment_mass is not None:
            for i, subset in enumerate(subsets_r2):
                subsets_r2[i] = [round(exact_mass - loss, 4)] + subset

        if debug:
            print("Second round (mass: {}) - Values: {} - Correct Sums: {}".format(exact_mass__0_0001,
                                                                                   len(mass_values_r2),
                                                                                   len(subsets_r2)))
            print("------------------------------------------------------")

        l = multiprocessing.Lock()
        with multiprocessing.Pool(initializer=multiprocess_init, initargs=(l,)) as pool:

            pool.map(partial(build_from_subsets, configs_iso=configs_iso, mc=mc, table_name=table_name, ppm=ppm,
                             fn_out=fn_out, debug=debug, path_db_k_graphs=path_db_k_graphs, path_pkls=path_pkls,
                             path_db=path_db),
                     subsets_r2)

    db.close()


def multiprocess_init(l):
    global lock
    lock = l


def gen_subs_table(db, heavy_atoms, max_valence, max_atoms_available):
    table_name = "subset_substructures"
    db.cursor.execute("DROP TABLE IF EXISTS %s" % table_name)
    db.cursor.execute("""CREATE TABLE {} AS
                             SELECT * FROM substructures WHERE
                                 heavy_atoms IN ({}) AND
                                 atoms_available <= {} AND
                                 valence <= {}
                      """.format(table_name,
                                 ",".join(map(str, heavy_atoms)),
                                 max_valence,
                                 max_atoms_available,))

    db.create_indexes(table=table_name, selection="gen_subs_table")

    return table_name


def build_from_subsets(ss2_grp, configs_iso, mc, table_name, fn_out,
                       path_db_k_graphs="../databases/k_graphs.sqlite", path_pkls="../databases/pkls",
                       path_db="../databases/substructures.sqlite", ppm=None, debug=False):
    """
    A stage of the :py:meth:`metaboverse.build_structures.build` workflow for generating molecules to a given mass
    from substructures. At this stage, mass subsets have been identified in the substructure database. Each of these
    groups are now filtered further by identifying masses that refer to valid subsets of molecules, before
    """

    out = open(fn_out, "a")

    db = SubstructureDb(path_db, path_pkls, path_db_k_graphs, False)
    list_ecs = combine_ecs(ss2_grp, db, table_name, "0_0001", ppm)

    if len(list_ecs) == 0:
        db.close()
        return

    iii = 0
    for l in itertools.product(*list_ecs):

        sum_ec = list(numpy.array(l).sum(axis=0))
        iii += 1

        if mc != sum_ec and debug:
            print("No match for elemental composition: {}".format(str(sum_ec)))

        elif mc == sum_ec:

            if debug:
                print("Match elemental composition: {}".format(str(sum_ec)))

            ll = db.select_sub_structures(l, table_name)

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

            if debug:
                print("## {} substructure combinations".format(len(list(itertools.product(*ll)))))

            smis = []
            for lll in itertools.product(*ll):

                if debug:
                    for record in lll:
                        print(record)
                    print("---------------")

                lll = sorted(lll, key=itemgetter('atoms_available', 'valence'))

                vA = ()
                for d in lll:
                     vA += (tuple(d["degree_atoms"].values()),)

                if str(vA) not in configs_iso:
                    if debug:
                        print("NO:", str(vA))
                        print("============")
                    continue
                else:
                    if debug:
                        print("YES:", str(vA))
                        print("============")

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
                        Chem.SanitizeMol(molOut)
                    except:
                        if debug:
                            print("Can't sanitize mol ISO: {}".format(iso_n))
                        continue

                    try:
                        smis.append("{}\t{}\n".format(Chem.MolToSmiles(molOut), str([item["smiles"] for item in lll])))
                    except RuntimeError:
                        if debug:
                            print("Bad bond type violation")

                    if debug:
                        print("## smi (result): {}".format(Chem.MolToSmiles(molOut)))

            with lock:
                for smi in smis:
                    out.write(smi)

    db.close()
