#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019-2020 Ralf Weber
#
# This file is part of MetaboBlend.
#
# MetaboBlend is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MetaboBlend is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MetaboBlend.  If not, see <https://www.gnu.org/licenses/>.
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


def find_path(l, dp, n, mass, max_subset_length, path=[]):
    """
    Recursive solution for backtracking through the dynamic programming boolean matrix. All possible subsets are found

    :param l: A list of masses from which to identify subsets.

    :param mass: The target mass of the sum of the substructures.

    :param dp: The dynamic programming boolean matrix.

    :param n: The size of l.

    :param max_subset_length: The maximum length of subsets to return. Allows the recursive backtracking algorithm to
        terminate early in many cases, significantly improving runtime.

    :param path: List for keeping track of the current subset.

    :return: Generates of lists containing the masses of valid subsets.
    """

    # base case - the path has generated a correct solution
    if mass == 0:
        yield sorted(path)
        return

    # stop running when we overshoot the mass
    elif mass < 0:
        return

    # can we sum up to the target value using the remaining masses? recursive call
    elif dp[n][mass]:
        yield from find_path(l, dp, n-1, mass, max_subset_length, path)

        if len(path) < max_subset_length:
            path.append(l[n-1])

            yield from find_path(l, dp, n-1, mass - l[n-1], max_subset_length, path)
            path.pop()


def subset_sum(l, mass, max_subset_length=3):
    """
    Dynamic programming implementation of subset sum. Note that, whilst this algorithm is pseudo-polynomial, the
    backtracking algorithm for obtaining all possible subsets has exponential complexity and so remains unsuitable
    for large input values.  This does, however, tend to perform a lot better than non-dp implementations, as we're
    no longer doing sums multiple times and we've cut down the operations performed during the exponential portion of
    the method.

    :param l: A list of masses from which to identify subsets.

    :param mass: The target mass of the sum of the substructures.

    :param max_subset_length: The maximum length of subsets to return. Allows the recursive backtracking algorithm to
        terminate early in many cases, significantly improving runtime.

    :return: Generates of lists containing the masses of valid subsets.
    """
    n = len(l)

    # initialise dynamic programming array
    dp = numpy.ndarray([n+1, mass+1], bool)

    # subsets can always equal 0
    for i in range(n+1):
        dp[i][0] = True

    # empty subsets do not have non-zero sums
    for i in range(mass):
        dp[0][i+1] = False

    # fill in the remaining boolean matrix
    for i in range(n):
        for j in range(mass+1):
            if j >= l[i]:
                dp[i+1][j] = dp[i][j] or dp[i][j-l[i]]
            else:
                dp[i+1][j] = dp[i][j]

    # backtrack through the matrix recursively to obtain all solutions
    return find_path(l, dp, n, mass, max_subset_length)


def combine_ecs(ss2_grp, db, table_name, accuracy, ppm=None):
    """
    A wrapper for :py:meth:`metaboblend.databases.select_ecs` that instead takes a group of subsets, as generated by
    the second stage of :py:meth:`metaboblend.build_structres.subset_sum` in
    :py:meth:`metaboblend.build_structres.build`.

    :param ss2_grp: A list containing the masses of substructures identified by subset_sum.

    :param db: The :py:meth:`metaboblend.databases.SubstructureDb` in which to search for elemental compositions.

    :param table_name: The name of the table containing substructures in which to search for elemental compositions.

    :param accuracy: To which decimal places of accuracy results are to be limited to.

            * **1** Integer level
            * **0_0001** Four decimal places

    :param ppm: The allowable error of the query (in parts per million). If unspecified, only exact matches are
        considered.

    :return: If there are no elemental compositions for any of the masses in the group, then an empty list is returned.
    """

    ecs = []

    for i in range(len(ss2_grp)):
        atoms = db.select_ecs(ss2_grp[i], table_name, accuracy, ppm=ppm)

        if len(atoms) == 0:
            return []

        ecs.append(atoms)

    return ecs


def reindex_atoms(records):
    """
    Parses the libs of groups of substructures that are to be combined; the lib is a dictionary containing details
    about the substructure, as generated by :py:meth:`metaboblend.databases.get_substructure`. Combines the
    molecules into a single :py:meth:`rdkit.Chem.Mol` object and obtains details on their bonding properties.

    :param records: Takes a list of lib dictionaries that contain details on each substructure to be combined.

    :return: Returns a tuple containing a :py:meth:`rdkit.Chem.CombineMols` object, that stores all the substructures
        as a single molecule, followed by information on the substructure bonding properties, including:

        * **atoms_available** A list of the indices of atoms that are available for bonding.

        * **atoms_to_remove** A list of the indices of dummy atoms that are to be removed in order to bond with other
        substructures.

        * **bond_types** A dictionary containing the indices of atoms that are available for bonding as keys and values
        detailing their bond types. See :py:meth:`metaboblend.build_structures.add_bonds`.
    """

    atoms_available, atoms_to_remove, bond_types = [], [], {}
    mol_comb = Chem.Mol()
    index_atoms, all_bond_types = [], {}
    c = 0

    for i, record in enumerate(records):
        idxs = []
        all_bond_types[i] = []
        for atom in record["mol"].GetAtoms():

            newIdx = atom.GetIdx() + c
            idxs.append(newIdx)

            if atom.GetIdx() in record["degree_atoms"]:
                atoms_available.append(newIdx)
            if atom.GetIdx() in record["dummies"]:
                atoms_to_remove.append(newIdx)
            if atom.GetIdx() in record["bond_types"]:
                bond_types[newIdx] = record["bond_types"][atom.GetIdx()]
                all_bond_types[i] += record["bond_types"][atom.GetIdx()]

        mol_comb = Chem.CombineMols(mol_comb, record["mol"])
        index_atoms.append(idxs)
        c = idxs[-1] + 1

    # check that bond types add up - removes some mismatched configurations
    bond_mismatch = False
    for i in range(len(records)):
        other_bonds = []
        for j in range(len(records)):
            if i != j:
                other_bonds += all_bond_types[j]

        for bond in all_bond_types[i]:
            if bond not in other_bonds:
                bond_mismatch = True

    return mol_comb, atoms_available, atoms_to_remove, bond_types, bond_mismatch


def add_bonds(mols, edges, atoms_available, bond_types, debug=False):
    """
    Takes a set of substructures and attempts to combine them together to generate a final structure. One of the last
    steps in the :py:meth:`metaboblend.build_structures.build` workflow.

    :param mols: A :py:meth:`rdkit.Chem.CombineMols` object, that stores all the substructures
        as a single molecule.

    :param edges: The edges to use in order to join the substructures together, obtained from the connectivity database
        (:py:meth:`metaboblend.databases.create_isomorphism_database`).

    :param atoms_available: A list of the indices of atoms that are available for bonding.

    :param bond_types:  The type of bonds to be formed by dummy atoms - see :py:meth:`Chem.rdchem.BondType`. Is a
        dictionary whose keys are atom indices and values are bond types, as follows:

            * **1.0** Single

            * **1.5** Aromatic

            * **2.0** Double

    :param debug: Debug print statements provide further information on how the function is generating the connectivity
        database.

        * **True** Print debug statements.

        * **False** Hide debug print statements.

    :return: If unsuccessful, returns None, else returns an :py:meth:`rdkit.Chem.EditableMol` object containing
        the substructures combined into a final single molecule.
    """

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
          path_db="../databases/substructures.sqlite", fragment_mass=None, ppm=None, debug=False, out_mode="w",
          processes=None, table_name=None, minimum_frequency=None):
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
            * **0_0001** Four decimal places

    :param max_n_substructures: The maximum number of substructures to be used for building molecules.

    :param path_db: The path to the SQLite 3 substructure database, as generated by
        :py:meth:`metaboblend.databases.SubstructureDb`.

    :param path_db_k_graphs: The path to the SQLite 3 connectivity database, as generated by
        :py:meth:`metaboblend.databases.create_isomorphism_database`.

    :param path_pkls: The path to the connectivity graphs described by the SQLite 3 connectivity database, as generated
        by :py:meth:`metaboblend.databases.create_isomorphism_database`.

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

    :param processes: How many worker processes to utilise; if left as None, :py:meth:`os.cpu_count` is used.

    :param table_name: If specified, the table specified within the substructure database will be used to generate
        molecules; else, a temporary substructure table will be created by the function.
    """

    db = SubstructureDb(path_db, path_pkls, path_db_k_graphs)

    if table_name is None:  # generate "temp" table containing only substructures in parameter space
        table_name = gen_subs_table(db, heavy_atoms, max_valence, max_atoms_available, round(exact_mass),
                                    minimum_frequency=minimum_frequency)

    if fragment_mass is None:  # standard build method
        exact_mass__1 = round(exact_mass)
        exact_mass__0_0001 = round(exact_mass, 4)

        tolerance = 0.001

        fragment_edges_only = False

    else:  # prescribed substructure build method
        loss = exact_mass - fragment_mass
        exact_mass__1 = round(loss)
        exact_mass__0_0001 = round(loss, 4)

        tolerance = (loss / 1000000) * ppm
        if tolerance < 0.001:
            tolerance = 0.001
        else:
            tolerance = round(tolerance, 4)

        max_n_substructures -= 1  # we find sets of mols that add up to the loss, not the precursor mass

        fragment_edges_only = True

    if os.name == "nt":  # multiprocessing freeze support on windows
        multiprocessing.freeze_support()

    # select groups of masses at low mass resolution
    mass_values = [m for m in db.select_mass_values("1", [], table_name) if m <= exact_mass__1]
    if len(mass_values) == 0:
        subsets = []
    else:
        subsets = list(subset_sum(mass_values, exact_mass__1, max_n_substructures))

    configs_iso = db.k_configs(fragment_edges_only)
    out = open(fn_out, out_mode)

    if debug:
        print("""First round (mass: {}) - Values: {} - Correct Sums: {}
              """.format(exact_mass__1, len(mass_values), len(subsets)))
        print("------------------------------------------------------")

    lls = []
    for ss_grp in subsets:
        if len(ss_grp) > max_n_substructures or len(ss_grp) == 0:
            continue

        # refine groups of masses to 4dp mass resolution
        mass_values_r2 = db.select_mass_values("0_0001", ss_grp, table_name)
        subsets_r2 = []

        # use combinations to get second group of masses instead of subset sum - subset sum is integer mass only
        for mass_combo in itertools.product(*mass_values_r2):
            if abs(sum(mass_combo) - exact_mass__0_0001) <= tolerance:
                subsets_r2.append(mass_combo)

        if len(subsets_r2) == 0:
            continue

        if fragment_mass is not None:  # add fragments mass to to loss group
            subsets_r2 = [subset + (round(exact_mass - loss, 4),) for subset in subsets_r2]

        if debug:
            print("""Second round (mass: {}) - Values: {} - Correct Sums: {}
                  """.format(exact_mass__0_0001, len(mass_values_r2), len(subsets_r2)))
            print("------------------------------------------------------")

        for ss_grp2 in subsets_r2:  # refines groups based on ecs and gets substructures from db (appends to lls)
            build_from_subsets(ss_grp2, mc=mc, table_name=table_name, ppm=ppm, debug=debug, db=db, lls=lls)

    with multiprocessing.Pool(processes=processes) as pool:  # send sets of substructures for building
        smi_list = pool.map(partial(lll_build, debug=debug, configs_iso=configs_iso), lls)

    if len(smi_list) != 0:
        out.writelines(smi_list)

    out.close()
    db.close()


def gen_subs_table(db, heavy_atoms, max_valence, max_atoms_available, max_mass, table_name="subset_substructures",
                   minimum_frequency=None):
    """
    Generate a temporary secondary substructure table restricted by a set of parameters. Generated as an initial step
    in :py:meth:`metaboblend.build_structures.build` in order to limit the processing overhead as a result of
    repeatedly querying the SQLite substructure database.

    :param max_mass: The maximum allowed mass of substructures in the temporary table; there is no point considering
        substructures with greater mass than the target mol.

    :param db: Connection to a :py:meth:`metaboblend.databases.SubstructureDb` from which to extract substructures.

    :param heavy_atoms: List of integers used to limit which substructures are transferred into the temporary table.

    :param max_valence: The maximum total valence (ie, the product of `atoms_available` and the degree of their bonds)
        to be included in the temporary table.

    :param max_atoms_available: The maximal atoms available of substructures to be included in the temporary table.

    :param table_name: Defaults to "subset_substructures", which is cleaned up upon database closure. The name of the
        table to be generated

    :return: The name of the temporary secondary substructure table.
    """

    db.cursor.execute("DROP TABLE IF EXISTS %s" % table_name)

    if minimum_frequency is None:
        freq_statement = ""
    else:
        freq_statement = """
                            AND smiles IN 
                                (SELECT smiles 
                                    FROM hmdbid_substructures 
                                    GROUP BY smiles 
                                    HAVING COUNT(*) >= {})
                            """.format(minimum_frequency,)

    db.cursor.execute("""CREATE TABLE {} AS
                             SELECT * FROM substructures WHERE
                                 heavy_atoms IN ({}) AND
                                 atoms_available <= {} AND
                                 valence <= {} AND
                                 exact_mass__1 < {}{}
                      """.format(table_name,
                                 ",".join(map(str, heavy_atoms)),
                                 max_atoms_available,
                                 max_valence,
                                 max_mass,
                                 freq_statement,))

    db.create_indexes(table=table_name, selection="gen_subs_table")

    return table_name


def build_from_subsets(ss2_grp, mc, table_name, db, lls=[], ppm=None, debug=False,):
    """
    A stage of the :py:meth:`metaboblend.build_structures.build` workflow for generating molecules to a given mass
    from substructures. At this stage, mass subsets have been identified in the substructure database. Each of these
    groups are now filtered further by identifying masses that refer to valid subsets of molecules, before they are
    built to generate new molecules.

    :param db: The substructure and connectivity database. Elemental compositions and substructures are retrieved from
        the database; this information is listed as "ll" and will be appended to the lls list provided as a parameter.

    :param lls: List of substructure combinations, as determined by this function. Note that this list will be appended
        to by the function, but the original items in the list will not be changed.

    :param ss2_grp: Group of masses that sum to the correct total mass, refer to substructures in the substructure
        database.

    :param mc: List of integers detailing the molecular composition of the target metabolite, in the format
        `[C, H, N, O, P, S]`.

    :param ppm: The allowable error of the query (in parts per million). Designed to be used in accordance with
        `fragment_mass`.

    :param debug: Debug print statements provide further information on how the function is generating the connectivity
        database.

        * **True** Print debug statements.

        * **False** Hide debug print statements.

    :param table_name: The name of the table within the substructure database from which to extract substructures. A
        prefiltered table based on the parameters specified in :py:meth:`metaboblend.build_structures.build`. See
        :py:meth:`metaboblend.build_structures.gen_subs_table`.
    """

    list_ecs = combine_ecs(ss2_grp, db, table_name, "0_0001", ppm)

    if len(list_ecs) == 0:
        return

    iii = 0
    for l in itertools.product(*list_ecs):

        sum_ec = list(numpy.array(l).sum(axis=0))
        iii += 1

        if mc != sum_ec:  # check each set of elemental compositions matches the target mol
            if debug:
                print("No match for elemental composition: {}".format(str(sum_ec)))

            continue

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

        lls += itertools.product(*ll)  # get the combinations of retrieved substructures


def lll_build(lll, configs_iso, debug):
    """
    Final stage for building molecules; takes a combination of substructures (lll) and builds them according to
    graphs in the substructure database. May be run in parallel.

    :param lll: Combinations of substructures for building mols.
        :param debug: Debug print statements provide further information on how the function is generating the connectivity
        database.

        * **True** Print debug statements.

        * **False** Hide debug print statements.

    :param configs_iso: Possible substructure combinations extracted from the connectivity database. A tuple containing
        tuples for each substructure; these tuples specify how many bonds each substructure can make.

    :return: List of smiles representing molecules generated (and the substructures used to generate them).
    """

    smis = ""

    if debug:
        for record in lll:
            print(record)
        print("---------------")

    lll = sorted(lll, key=itemgetter('atoms_available', 'valence'))

    vA = ()
    for d in lll:
        vA += (tuple(d["degree_atoms"].values()),)  # obtain valence configuration of the set of substructures

    if str(vA) not in configs_iso:  # check mols "fit" together according to the connectivity database
        if debug:
            print("NO:", str(vA))
            print("============")
        return ""

    else:
        if debug:
            print("YES:", str(vA))
            print("============")

    mol_comb, atoms_available, atoms_to_remove, bond_types, bond_mismatch = reindex_atoms(lll)

    if bond_mismatch:
        return ""  # check that bond types are compatible (imperfect check)

    if debug:
        print("## Mols (in memory):", mol_comb)
        print("## Atoms Available (indexes):", atoms_available)
        print("## Atoms to remove (dummies):", atoms_to_remove)
        print("## Type of bonds to form:", bond_types)

    iso_n = 0
    for edges in configs_iso[str(vA)]:  # build mols for each graph in connectivity db
        iso_n += 1
        if debug:
            print("## ISO {}".format(iso_n))

        if debug:
            print("1: Add bonds")

        mol_e = add_bonds(mol_comb, edges, atoms_available, bond_types)  # add bonds between substructures

        if mol_e is None:
            continue
        if debug:
            print("2: Add bonds")

        atoms_to_remove.sort(reverse=True)
        [mol_e.RemoveAtom(a) for a in atoms_to_remove]  # clean up dummy atoms

        molOut = mol_e.GetMol()  # generate the final (non-editable) mol

        try:
            Chem.SanitizeMol(molOut)  # clean the mol - ensure it is valid & canonical
        except:
            if debug:
                print("Can't sanitize mol ISO: {}".format(iso_n))
            continue

        try:  # append the canonical smiles of the final structure
            smis += "{}\n".format(Chem.MolToSmiles(molOut))
        except RuntimeError:
            if debug:
                print("Bad bond type violation")

        if debug:
            print("## smi (result): {}".format(Chem.MolToSmiles(molOut)))

    return smis
