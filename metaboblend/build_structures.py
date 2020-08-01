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


def find_path(mass_list, sum_matrix, n, mass, max_subset_length, path=[]):
    """
    Recursive solution for backtracking through the dynamic programming boolean matrix. All possible subsets are found

    :param mass_list: A list of masses from which to identify subsets.

    :param mass: The target mass of the sum of the substructures.

    :param sum_matrix: The dynamic programming boolean matrix.

    :param n: The size of mass_list.

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
    elif sum_matrix[n][mass]:
        yield from find_path(mass_list, sum_matrix, n - 1, mass, max_subset_length, path)

        if len(path) < max_subset_length:
            path.append(mass_list[n-1])

            yield from find_path(mass_list, sum_matrix, n - 1, mass - mass_list[n - 1], max_subset_length, path)
            path.pop()


def subset_sum(mass_list, mass, max_subset_length=3):
    """
    Dynamic programming implementation of subset sum. Note that, whilst this algorithm is pseudo-polynomial, the
    backtracking algorithm for obtaining all possible subsets has exponential complexity and so remains unsuitable
    for large input values.  This does, however, tend to perform a lot better than non-sum_matrix implementations, as we're
    no longer doing sums multiple times and we've cut down the operations performed during the exponential portion of
    the method.

    :param mass_list: A list of masses from which to identify subsets.

    :param mass: The target mass of the sum of the substructures.

    :param max_subset_length: The maximum length of subsets to return. Allows the recursive backtracking algorithm to
        terminate early in many cases, significantly improving runtime.

    :return: Generates of lists containing the masses of valid subsets.
    """

    n = len(mass_list)

    # initialise dynamic programming array
    sum_matrix = numpy.ndarray([n + 1, mass + 1], bool)

    # subsets can always equal 0
    for i in range(n+1):
        sum_matrix[i][0] = True

    # empty subsets do not have non-zero sums
    for i in range(mass):
        sum_matrix[0][i + 1] = False

    # fill in the remaining boolean matrix
    for i in range(n):
        for j in range(mass+1):
            if j >= mass_list[i]:
                sum_matrix[i + 1][j] = sum_matrix[i][j] or sum_matrix[i][j - mass_list[i]]
            else:
                sum_matrix[i + 1][j] = sum_matrix[i][j]

    # backtrack through the matrix recursively to obtain all solutions
    return find_path(mass_list, sum_matrix, n, mass, max_subset_length)


def combine_ecs(precise_mass_grp, db, table_name, accuracy):
    """
    A wrapper for :py:meth:`metaboblend.databases.select_ecs` that instead takes a group of subsets, as generated by
    the second stage of :py:meth:`metaboblend.build_structres.subset_sum` in
    :py:meth:`metaboblend.build_structres.build`.

    :param precise_mass_grp: A list containing the masses of substructures identified by subset_sum.

    :param db: The :py:meth:`metaboblend.databases.SubstructureDb` in which to search for elemental compositions.

    :param table_name: The name of the table containing substructures in which to search for elemental compositions.

    :param accuracy: To which decimal places of accuracy results are to be limited to.

            * **1** Integer level
            * **0_0001** Four decimal places

    :return: If there are no elemental compositions for any of the masses in the group, then an empty list is returned.
    """

    ecs = []

    for i in range(len(precise_mass_grp)):
        atoms = db.select_ecs(precise_mass_grp[i], table_name, accuracy)

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

            new_idx = atom.GetIdx() + c
            idxs.append(new_idx)

            if atom.GetIdx() in record["degree_atoms"]:
                atoms_available.append(new_idx)
            if atom.GetIdx() in record["dummies"]:
                atoms_to_remove.append(new_idx)
            if atom.GetIdx() in record["bond_types"]:
                bond_types[new_idx] = record["bond_types"][atom.GetIdx()]
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


def add_bonds(mols, edges, atoms_available, bond_types):
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

    :return: If unsuccessful, returns None, else returns an :py:meth:`rdkit.Chem.EditableMol` object containing
        the substructures combined into a final single molecule.
    """

    rdkit_bond_types = {1: Chem.rdchem.BondType.SINGLE,
                        1.5: Chem.rdchem.BondType.AROMATIC,
                        2: Chem.rdchem.BondType.DOUBLE}

    g = nx.Graph()
    g.add_edges_from(edges)

    g = nx.relabel_nodes(g, dict(zip(sorted(g.nodes()), atoms_available)))

    mol_edit = Chem.EditableMol(mols)
    for edge in g.edges():

        if edge[0] in bond_types:
            bt_start = copy.copy(bond_types[edge[0]])
        else:
            return None  # nested dummy

        if edge[1] in bond_types:
            bt_end = copy.copy(bond_types[edge[1]])
        else:
            return None  # nested dummy

        bond_matches = list(set(bt_start).intersection(bt_end))

        if len(bond_matches) == 0:
            return None

        else:
            bt_start.remove(bond_matches[0])
            bt_end.remove(bond_matches[0])

        try:
            mol_edit.AddBond(edge[0], edge[1], rdkit_bond_types[bond_matches[0]])
        except KeyError:
            return None  # unknown bond type

    return mol_edit


def annotate_msn(mc, exact_mass, fragment_masses, smi_out_dir=None, heavy_atoms=range(0, 10), max_valence=6,
                 max_atoms_available=2, max_n_substructures=3, path_connectivity_db="../databases/k_graphs.sqlite",
                 path_substructure_db="../databases/substructures.sqlite", ppm=5, processes=None,
                 write_fragment_smis=False, return_smi_dict=False, minimum_frequency=None, hydrogenation_allowance=2):
    """
    :param mc: List of integers detailing the molecular composition of the target metabolite, in the format
        [C, H, N, O, P, S]. Can also be a list of elemental composition lists for generating multiple structures; in
        this case, exact_mass and fragment_masses must also be a list of the same length.

    :param exact_mass: The exact mass of the target metabolite, or list of masses.

    :param fragment_masses: A list of the neutral masses generated at a higher order MSn level for predict the structure
        for the supplied mc/exact_mass (MSn-1). Can also be a list of lists for generating multiple structures.

    :param smi_out_dir: The path of the file to which unique smile strings should be written representing the final
        structures generated. If None, no file is written.

    :param heavy_atoms: A list of integers containing the heavy atom counts to consider in substructures to be used
        to build final structures.

    :param max_valence: The maximal total bond orders of substructures to be considered to build final structures
        (ie, the product of `atoms_available` and the degree of their bonds).

    :param max_atoms_available: The maximal atoms available of substructures to be considered for building molecules.

    :param max_n_substructures: The maximum number of substructures to be used for building molecules.

    :param path_substructure_db: The path to the SQLite 3 substructure database, as generated by
        :py:meth:`metaboblend.databases.SubstructureDb`.

    :param path_connectivity_db: The path to the SQLite 3 connectivity database, as generated by
        :py:meth:`metaboblend.databases.create_isomorphism_database`.

    :param processes: How many worker processes to utilise; if left as None, :py:meth:`os.cpu_count` is used.

    :param minimum_frequency: The minimum frequency of substructures in table_name; e.g. substructures have a frequency
        of 1 if they are unique.

    :param hydrogenation_allowance:

    :param ppm: The allowable error of the query (in parts per million). Designed to be used in accordance with
        `fragment_mass`.

    :param write_fragment_smis: Whether to write smiles to a file for each fragment.

    :param return_smi_dict: Whether to return a dict of smiles.
    """

    if isinstance(mc[1], int) and isinstance(exact_mass, int):  # single input
        mc, exact_mass = [mc], [exact_mass]
        multi_input = False
    elif isinstance(mc[1], list) and isinstance(exact_mass, list):  # multiple input
        multi_input = True
    else:
        raise ValueError("either pass a single input to mc and exact_mass, or lists of the same length")

    # prepare temporary table here - will only be generated once in case of multiple input
    db = SubstructureDb(path_substructure_db, path_connectivity_db)
    table_name = gen_subs_table(db, heavy_atoms, max_valence, max_atoms_available, round(max(exact_mass)),
                                minimum_frequency=minimum_frequency)
    db.close()

    for i, curr_mc, curr_exact_mass, curr_fragment_masses in zip(range(len(mc)), mc, exact_mass, fragment_masses):
        if smi_out_dir is not None and write_fragment_smis:
            smi_out_subdir = os.path.join(smi_out_dir, str(i) + "_" + str(round(curr_exact_mass)))
            os.mkdir(smi_out_subdir)
        else:
            smi_out_subdir = None

        smi_dict = build_msn(mc=curr_mc, exact_mass=curr_exact_mass, fragment_masses=curr_fragment_masses,
                             heavy_atoms=heavy_atoms, max_valence=max_valence, max_atoms_available=max_atoms_available,
                             max_n_substructures=max_n_substructures, smi_out_dir=smi_out_subdir,
                             path_connectivity_db=path_connectivity_db, path_substructure_db=path_substructure_db,
                             processes=processes, minimum_frequency=minimum_frequency, return_smi_dict=return_smi_dict,
                             hydrogenation_allowance=hydrogenation_allowance, ppm=ppm, table_name=table_name,
                             write_fragment_smis=write_fragment_smis)

        if return_smi_dict and multi_input:
            yield smi_dict

    if return_smi_dict and not multi_input:
        return smi_dict


def build_msn(mc, exact_mass, fragment_masses, heavy_atoms, max_valence, max_atoms_available, max_n_substructures,
              smi_out_dir, path_connectivity_db, path_substructure_db, processes, minimum_frequency, return_smi_dict,
              hydrogenation_allowance, ppm, table_name, write_fragment_smis):
    """
    Additional considerations for building from MS/MS data.

    :param mc: List of integers detailing the molecular composition of the target metabolite, in the format
        [C, H, N, O, P, S]. Can also be a list of elemental composition lists for generating multiple structures; in
        this case, exact_mass and fragment_masses must also be a list of the same length.

    :param exact_mass: The exact mass of the target metabolite, or list of masses.

    :param fragment_masses: A list of the neutral masses generated at a higher order MSn level for predict the structure
        for the supplied mc/exact_mass (MSn-1). Can also be a list of lists for generating multiple structures.

    :param smi_out_dir: The path of the file to which unique smile strings should be written representing the final
        structures generated. If None, no file is written.

    :param heavy_atoms: A list of integers containing the heavy atom counts to consider in substructures to be used
        to build final structures.

    :param max_valence: The maximal total bond orders of substructures to be considered to build final structures
        (ie, the product of `atoms_available` and the degree of their bonds).

    :param max_atoms_available: The maximal atoms available of substructures to be considered for building molecules.

    :param max_n_substructures: The maximum number of substructures to be used for building molecules.

    :param path_substructure_db: The path to the SQLite 3 substructure database, as generated by
        :py:meth:`metaboblend.databases.SubstructureDb`.

    :param path_connectivity_db: The path to the SQLite 3 connectivity database, as generated by
        :py:meth:`metaboblend.databases.create_isomorphism_database`.

    :param processes: How many worker processes to utilise; if left as None, :py:meth:`os.cpu_count` is used.

    :param minimum_frequency: The minimum frequency of substructures in table_name; e.g. substructures have a frequency
        of 1 if they are unique.

    :param hydrogenation_allowance: Allowable number of hydrogen masses that each fragment can deviate from the "true"
        substructure.

    :param ppm: The allowable error of the query (in parts per million). Designed to be used in accordance with
        `fragment_mass`.

    :param table_name: If specified, the table specified within the substructure database will be used to generate
        molecules; else, a temporary substructure table will be created by the function.

    :param write_fragment_smis: Whether to write smiles to a file for each fragment.

    :param return_smi_dict: Whether to return a dict of smiles.
    """
    structure_frequency = {}  # map smiles to how many separate masses generated them
    fragment_masses.sort(reverse=True)

    if table_name is None:
        db = SubstructureDb(path_substructure_db, path_connectivity_db)
        table_name = gen_subs_table(db, heavy_atoms, max_valence, max_atoms_available, round(max(exact_mass)),
                                    minimum_frequency=minimum_frequency)
        db.close()

    for fragment_mass in fragment_masses:
        if smi_out_dir is not None and write_fragment_smis:
            smi_out = os.path.join(smi_out_dir, str(round(fragment_mass, 4)) + ".smi")
            open(smi_out, "w").close()
        else:
            smi_out = None

        fragment_smis = set()

        for j in range(0 - hydrogenation_allowance, hydrogenation_allowance + 1):
            hydrogenated_fragment_mass = fragment_mass + (j * 1.007825)  # consider re-arrangements
            fragment_smis.update(build(mc=mc, exact_mass=exact_mass, smi_out=smi_out, heavy_atoms=heavy_atoms,
                                       max_valence=max_valence, max_atoms_available=max_atoms_available,
                                       max_n_substructures=max_n_substructures,
                                       path_connectivity_db=path_connectivity_db,
                                       path_substructure_db=path_substructure_db,
                                       fragment_mass=hydrogenated_fragment_mass,
                                       ppm=ppm, out_mode="a", table_name=table_name,
                                       processes=processes, minimum_frequency=None))

        for smi in fragment_smis:
            structure_frequency[smi] = structure_frequency.get(smi, 0) + 1

    if smi_out_dir is not None:
        with open(os.path.join(smi_out_dir, "structure_frequency.csv"), "w") as freq_out:
            freq_out.writelines([k + "," + str(i) for k, i in zip(structure_frequency.keys(), structure_frequency.values())])

    if return_smi_dict:
        return structure_frequency


def generate_structures(mc, exact_mass, heavy_atoms=range(2,9), max_valence=6, max_atoms_available=2,
                        max_n_substructures=3, smi_out=None, path_connectivity_db="../databases/k_graphs.sqlite",
                        path_substructure_db="../databases/substructures.sqlite", prescribed_substructure_mass=None,
                        processes=None, minimum_frequency=None, return_smi_list=True):
    """
    Workflow for generating molecules of a given mass using substructures and connectivity graphs. Can optionally
    take a "prescribed" fragment mass to further filter results. Final structures are returned as a list and/or
    written in text format.

    :param mc: List of integers detailing the molecular composition of the target metabolite, in the format
        [C, H, N, O, P, S]. Can also be a list of elemental composition lists for generating multiple structures; in
        this case, exact_mass must also be a list of the same length.

    :param exact_mass: The exact mass of the target metabolite, or list of masses.

    :param smi_out: The path of the file to which unique smile strings should be written representing the final
        structures generated. If None, no file is written.

    :param heavy_atoms: A list of integers containing the heavy atom counts to consider in substructures to be used
        to build final structures.

    :param max_valence: The maximal total bond orders of substructures to be considered to build final structures
        (ie, the product of `atoms_available` and the degree of their bonds).

    :param max_atoms_available: The maximal atoms available of substructures to be considered for building molecules.

    :param max_n_substructures: The maximum number of substructures to be used for building molecules.

    :param path_substructure_db: The path to the SQLite 3 substructure database, as generated by
        :py:meth:`metaboblend.databases.SubstructureDb`.

    :param path_connectivity_db: The path to the SQLite 3 connectivity database, as generated by
        :py:meth:`metaboblend.databases.create_isomorphism_database`.

    :param processes: How many worker processes to utilise; if left as None, :py:meth:`os.cpu_count` is used.

    :param prescribed_substructure_mass: Mass of a prescribed substructure to subset results by, similar to OMG's
        system.

    :param minimum_frequency: The minimum frequency of substructures in table_name; e.g. substructures have a frequency
        of 1 if they are unique.

    :param return_smi_list: Whether to yield a list of smiles for each input.

    :return: Returns a list of unique smiles.
    """

    if isinstance(mc[1], int) and isinstance(exact_mass, int):  # single input
        mc, exact_mass = [mc], [exact_mass]
        multi_input = False
    elif isinstance(mc[1], list) and isinstance(exact_mass, list):  # multiple input
        multi_input = True
    else:
        raise ValueError("either pass a single input to mc and exact_mass, or lists of the same length")

    # prepare temporary table here - will only be generated once in case of multiple input
    db = SubstructureDb(path_substructure_db, path_connectivity_db)
    table_name = gen_subs_table(db, heavy_atoms, max_valence, max_atoms_available, round(max(exact_mass)),
                                minimum_frequency=minimum_frequency)
    db.close()

    for curr_mc, curr_exact_mass in zip(mc, exact_mass):
        smi_list = build(mc=curr_mc, exact_mass=curr_exact_mass, heavy_atoms=heavy_atoms, max_valence=max_valence,
                         max_atoms_available=max_atoms_available, max_n_substructures=max_n_substructures,
                         smi_out=smi_out, path_connectivity_db=path_connectivity_db,
                         path_substructure_db=path_substructure_db, fragment_mass=prescribed_substructure_mass,
                         ppm=None, out_mode="w", table_name=table_name, minimum_frequency=minimum_frequency,
                         processes=processes)

        if return_smi_list and multi_input:
            yield smi_list

    if return_smi_list and not multi_input:
        return smi_list


def build(mc, exact_mass, heavy_atoms, max_valence, max_atoms_available, max_n_substructures, smi_out,
          path_connectivity_db, path_substructure_db, fragment_mass, ppm, out_mode, processes, table_name,
          minimum_frequency):
    """
    Workflow for generating molecules of a given mass using substructures and connectivity graphs. Can optionally
    take a "prescribed" fragment mass to further filter results; this can be used to incorporate MSn data. Final
    molecules are written to the specified file in smiles format, which the substructures utilised to generate them.

    :param mc: List of integers detailing the molecular composition of the target metabolite, in the format
        [C, H, N, O, P, S].

    :param exact_mass: The exact mass of the target metabolite.

    :param smi_out: The path of the file to which unique smile strings should be written representing the final
        structures generated. If None, no file is written.

    :param heavy_atoms: A list of integers containing the heavy atom counts to consider in substructures to be used
        to build final structures.

    :param max_valence: The maximal total bond orders of substructures to be considered to build final structures
        (ie, the product of `atoms_available` and the degree of their bonds).

    :param max_atoms_available: The maximal atoms available of substructures to be considered for building molecules.

    :param max_n_substructures: The maximum number of substructures to be used for building molecules.

    :param path_substructure_db: The path to the SQLite 3 substructure database, as generated by
        :py:meth:`metaboblend.databases.SubstructureDb`.

    :param path_connectivity_db: The path to the SQLite 3 connectivity database, as generated by
        :py:meth:`metaboblend.databases.create_isomorphism_database`.

    :param fragment_mass: A mass by which to filter results; if not provided, all possible structures will be generated.

    :param ppm: The allowable error of the query (in parts per million). Designed to be used in accordance with
        `fragment_mass`.

    :param out_mode: The mode in which to write to the output file.

        * **"w"** Create a new file, overwriting any existing file.

        * **"w"** Append results to an existing file.

    :param processes: How many worker processes to utilise; if left as None, :py:meth:`os.cpu_count` is used.

    :param table_name: If specified, the table specified within the substructure database will be used to generate
        molecules; else, a temporary substructure table will be created by the function.

    :param minimum_frequency: The minimum frequency of substructures in table_name; e.g. substructures have a frequency
        of 1 if they are unique.
    """

    db = SubstructureDb(path_substructure_db, path_connectivity_db)

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
    if smi_out is not None:
        out = open(smi_out, out_mode)

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

        for ss_grp2 in subsets_r2:  # refines groups based on ecs and gets substructures from db (appends to lls)
            lls += build_from_subsets(ss_grp2, mc=mc, table_name=table_name, db=db)

    with multiprocessing.Pool(processes=processes) as pool:  # send sets of substructures for building
        smi_lists = pool.map(partial(lll_build, configs_iso=configs_iso), lls)

    smis = set([val for sublist in smi_lists for val in sublist])

    if smi_out is not None:
        if len(smis) != 0:
            out.writelines(smis)

        out.close()

    db.close()

    return smis


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

    :param minimum_frequency: The minimum frequency of substructures in table_name; e.g. substructures have a frequency
        of 1 if they are unique.
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


def build_from_subsets(ss2_grp, mc, table_name, db):
    """
    A stage of the :py:meth:`metaboblend.build_structures.build` workflow for generating molecules to a given mass
    from substructures. At this stage, mass subsets have been identified in the substructure database. Each of these
    groups are now filtered further by identifying masses that refer to valid subsets of molecules, before they are
    built to generate new molecules.

    :param db: The substructure and connectivity database. Elemental compositions and substructures are retrieved from
        the database; this information is listed as "ll" and will be appended to the lls list provided as a parameter.

    :param ss2_grp: Group of masses that sum to the correct total mass, refer to substructures in the substructure
        database.

    :param mc: List of integers detailing the molecular composition of the target metabolite, in the format
        `[C, H, N, O, P, S]`.

    :param table_name: The name of the table within the substructure database from which to extract substructures. A
        prefiltered table based on the parameters specified in :py:meth:`metaboblend.build_structures.build`. See
        :py:meth:`metaboblend.build_structures.gen_subs_table`.
    """

    lls = []

    list_ecs = combine_ecs(ss2_grp, db, table_name, "0_0001")

    if len(list_ecs) == 0:
        return []

    iii = 0
    for l in itertools.product(*list_ecs):

        sum_ec = list(numpy.array(l).sum(axis=0))
        iii += 1

        if mc != sum_ec:  # check each set of elemental compositions matches the target mol
            continue

        ll = db.select_sub_structures(l, table_name)

        if len(ll) == 0:
            continue

        lls.append(ll)  # get the combinations of retrieved substructures

    return lls


def lll_build(ll, configs_iso):
    """
    Final stage for building molecules; takes a combination of substructures (lll) and builds them according to
    graphs in the substructure database. May be run in parallel.

    :param ll: Combinations of substructures for building mols.

    :param configs_iso: Possible substructure combinations extracted from the connectivity database. A tuple containing
        tuples for each substructure; these tuples specify how many bonds each substructure can make.

    :return: List of smiles representing molecules generated (and the substructures used to generate them).
    """

    smis = []

    for lll in itertools.product(*ll):
        lll = sorted(lll, key=itemgetter('atoms_available', 'valence'))

        v_a = ()
        for d in lll:
            v_a += (tuple(d["degree_atoms"].values()),)  # obtain valence configuration of the set of substructures

        if str(v_a) not in configs_iso:  # check mols "fit" together according to the connectivity database
            continue

        mol_comb, atoms_available, atoms_to_remove, bond_types, bond_mismatch = reindex_atoms(lll)

        if bond_mismatch:
            continue  # check that bond types are compatible (imperfect check)

        for edges in configs_iso[str(v_a)]:  # build mols for each graph in connectivity db
            mol_e = add_bonds(mol_comb, edges, atoms_available, bond_types)  # add bonds between substructures

            if mol_e is None:
                continue

            atoms_to_remove.sort(reverse=True)
            [mol_e.RemoveAtom(a) for a in atoms_to_remove]  # clean up dummy atoms

            mol_out = mol_e.GetMol()  # generate the final (non-editable) mol

            try:
                Chem.SanitizeMol(mol_out)  # clean the mol - ensure it is valid & canonical
            except:
                continue

            try:  # append the canonical smiles of the final structure
                smis.append(Chem.MolToSmiles(mol_out))
            except RuntimeError:
                continue  # bad bond type violation

    return smis
