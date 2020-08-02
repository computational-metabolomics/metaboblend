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


def annotate_msn(ms_data, smi_out_dir=None, heavy_atoms=range(0, 10), max_valence=6, max_atoms_available=2,
                 max_n_substructures=3, path_connectivity_db="../databases/k_graphs.sqlite",
                 path_substructure_db="../databases/substructures.sqlite", ppm=5, processes=None,
                 write_fragment_smis=False, yield_smi_dict=True, minimum_frequency=None, hydrogenation_allowance=2,
                 add_small_substructures=True):
    """
    :param ms_data: Dictionary in the form ms_data[id] =
        {mc: [C, H, N, O, P, S], exact_mass: exact_mass, prescribed_masses=[]}. id represents a unique identifier for
        a given test, mc is a list of integers referring to molecular composition of the structure of interest,
        exact_mass is the mass of this structure to >=4d.p. and prescribed_masses are neutral fragment masses generated
        by this structure used to inform candidate scoring.

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

    :param hydrogenation_allowance: In order to represent re-arrangement events (the movement of hydrogens), in addition
        to attempting to build from substructures in prescribed_masses, we also attempt to build from prescribed_masses
        +- hydrogenation_allowance. E.g. if prescribed_masses = [2] then, to find candidate fragment substructures,
        we use as query the masses [2 - 1.007825, 2, 2 + 1.007825].

    :param ppm: The allowable error of the query (in parts per million).

    :param write_fragment_smis: Whether to write smiles to a file for each fragment.

    :param yield_smi_dict: Whether to return a dict of smiles.

    :param add_small_substructures: Whether to add a curated list of small substructures (single atom) to the filtered
        database prior to annotation.
    """

    if ppm is None:
        ppm = 0

    db = SubstructureDb(path_substructure_db, path_connectivity_db)

    # prepare temporary table here - will only be generated once in case of multiple input
    table_name = gen_subs_table(
        db=db,
        heavy_atoms=heavy_atoms,
        max_valence=max_valence,
        max_atoms_available=max_atoms_available,
        minimum_frequency=minimum_frequency,
        max_mass=round(max([ms_data[ms_id]["exact_mass"] for ms_id in ms_data.keys()]))
    )

    if add_small_substructures:
        db.add_small_substructures(table_name)

    for i, ms_id in enumerate(ms_data.keys()):
        if smi_out_dir is not None and write_fragment_smis:
            smi_out_path_subdir = os.path.join(smi_out_dir, str(i) + "_" + str(round(ms_data[ms_id]["exact_mass"])))
            os.mkdir(smi_out_path_subdir)
        else:
            smi_out_path_subdir = None

        smi_dict = build_msn(
            mc=ms_data[ms_id]["mc"],
            exact_mass=ms_data[ms_id]["exact_mass"],
            prescribed_masses=ms_data[ms_id]["prescribed_masses"],
            max_n_substructures=max_n_substructures,
            smi_out_dir=smi_out_path_subdir,
            path_connectivity_db=path_connectivity_db,
            path_substructure_db=path_substructure_db,
            processes=processes,
            return_smi_dict=yield_smi_dict,
            hydrogenation_allowance=hydrogenation_allowance,
            ppm=ppm,
            table_name=table_name,
            write_fragment_smis=write_fragment_smis
        )

        if yield_smi_dict:
            yield smi_dict

    db.close()


def build_msn(mc, exact_mass, prescribed_masses, max_n_substructures, smi_out_dir, path_connectivity_db,
              path_substructure_db, processes, return_smi_dict, hydrogenation_allowance, ppm, table_name,
              write_fragment_smis):
    """
    Additional considerations for building from MS/MS data.

    :param mc: List of integers detailing the molecular composition of the target metabolite, in the format
        [C, H, N, O, P, S]. Can also be a list of elemental composition lists for generating multiple structures; in
        this case, exact_mass and prescribed_masses must also be a list of the same length.

    :param exact_mass: The exact mass of the target metabolite, or list of masses.

    :param prescribed_masses: A list of the neutral masses generated at a higher order MSn level for predict the structure
        for the supplied mc/exact_mass (MSn-1). Can also be a list of lists for generating multiple structures.

    :param smi_out_dir: The path of the file to which unique smile strings should be written representing the final
        structures generated. If None, no file is written.

    :param max_n_substructures: The maximum number of substructures to be used for building molecules.

    :param path_substructure_db: The path to the SQLite 3 substructure database, as generated by
        :py:meth:`metaboblend.databases.SubstructureDb`.

    :param path_connectivity_db: The path to the SQLite 3 connectivity database, as generated by
        :py:meth:`metaboblend.databases.create_isomorphism_database`.

    :param processes: How many worker processes to utilise; if left as None, :py:meth:`os.cpu_count` is used.

    :param hydrogenation_allowance: In order to represent re-arrangement events (the movement of hydrogens), in addition
        to attempting to build from substructures in prescribed_masses, we also attempt to build from prescribed_masses
        +- hydrogenation_allowance. E.g. if prescribed_masses = [2] and hydrogenation_allowance = 1 then, to find
        candidate fragment substructures, we use as query the masses [2 - 1.007825, 2, 2 + 1.007825].

    :param ppm: The allowable error of the query (in parts per million). Designed to be used in accordance with
        `prescribed_mass`.

    :param table_name: If specified, the table specified within the substructure database will be used to generate
        molecules; else, a temporary substructure table will be created by the function.

    :param write_fragment_smis: Whether to write smiles to a file for each fragment.

    :param return_smi_dict: Whether to return a dict of smiles.
    """

    structure_frequency = {}  # map smiles to how many separate masses generated them
    prescribed_masses.sort(reverse=True)

    for prescribed_mass in prescribed_masses:
        if smi_out_dir is not None and write_fragment_smis:
            smi_out_path = os.path.join(smi_out_dir, str(round(prescribed_mass, 4)) + ".smi")
            open(smi_out_path, "w").close()
        else:
            smi_out_path = None

        fragment_smis = set()

        for j in range(0 - hydrogenation_allowance, hydrogenation_allowance + 1):
            hydrogenated_prescribed_mass = prescribed_mass + (j * 1.007825)  # consider re-arrangements

            fragment_smis.update(build(
                mc=mc,
                exact_mass=exact_mass,
                smi_out_path=smi_out_path,
                max_n_substructures=max_n_substructures,
                path_connectivity_db=path_connectivity_db,
                path_substructure_db=path_substructure_db,
                prescribed_mass=hydrogenated_prescribed_mass,
                ppm=ppm,
                out_mode="a",
                table_name=table_name,
                processes=processes,
                clean=False
            ))

        for smi in fragment_smis:
            structure_frequency[smi] = structure_frequency.get(smi, 0) + 1

    # write structure_frequency dict as csv
    if smi_out_dir is not None:
        with open(os.path.join(smi_out_dir, "structure_frequency.csv"), "w") as freq_out:
            freq_out.writelines([k + "," + str(i) + "\n" for k, i in zip(structure_frequency.keys(), structure_frequency.values())])

    if return_smi_dict:
        return structure_frequency


def generate_structures(ms_data, heavy_atoms=range(2, 9), max_valence=6, max_atoms_available=2,
                        max_n_substructures=3, smi_out_dir=None, path_connectivity_db="../databases/k_graphs.sqlite",
                        path_substructure_db="../databases/substructures.sqlite", processes=None,
                        minimum_frequency=None, yield_smi_list=True):
    """
    Workflow for generating molecules of a given mass using substructures and connectivity graphs. Can optionally
    take a "prescribed" fragment mass to further filter results. Final structures are returned as a list and/or
    written in text format.

    :param ms_data: Dictionary in the form ms_data[id] = 
        {mc: [C, H, N, O, P, S], exact_mass: structure mass, prescribed_masses: substructure mass}.

    :param smi_out_dir: The directory to which unique smile strings should be written representing the final
        structures generated. If None, no files are written.

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

    :param yield_smi_list: Whether to yield a list of smiles for each input.

    :return: Returns a list of unique smiles.
    """

    db = SubstructureDb(path_substructure_db, path_connectivity_db)
    smi_out_path = None

    # prepare temporary table here - will only be generated once in case of multiple input
    table_name = gen_subs_table(
        db=db,
        heavy_atoms=heavy_atoms,
        max_valence=max_valence,
        max_atoms_available=max_atoms_available,
        minimum_frequency=minimum_frequency,
        max_mass=round(max([ms_data[ms_id]["exact_mass"] for ms_id in ms_data.keys()]))
    )

    for ms_id in ms_data.keys():

        ppm = None
        try:
            if ms_data[ms_id]["prescribed_masses"] is not None:
                ppm = 0
        except KeyError:
            ms_data[ms_id]["prescribed_masses"] = None

        if smi_out_dir is not None:
            smi_out_path = os.path.join(smi_out_dir, ms_id + ".smi")

        smi_list = build(
            mc=ms_data[ms_id]["mc"],
            exact_mass=ms_data[ms_id]["exact_mass"],
            max_n_substructures=max_n_substructures,
            smi_out_path=smi_out_path,
            path_connectivity_db=path_connectivity_db,
            path_substructure_db=path_substructure_db,
            prescribed_mass=ms_data[ms_id]["prescribed_masses"],
            ppm=ppm,
            out_mode="w",
            table_name=table_name,
            processes=processes,
            clean=False
        )

        if yield_smi_list:
            yield smi_list

    db.close()


def build(mc, exact_mass, max_n_substructures, smi_out_path, path_connectivity_db, path_substructure_db,
          prescribed_mass, ppm, out_mode, processes, table_name, clean):
    """
    Workflow for generating molecules of a given mass using substructures and connectivity graphs. Can optionally
    take a "prescribed" fragment mass to further filter results; this can be used to incorporate MSn data. Final
    molecules are written to the specified file in smiles format, which the substructures utilised to generate them.

    :param mc: List of integers detailing the molecular composition of the target metabolite, in the format
        [C, H, N, O, P, S].

    :param exact_mass: The exact mass of the target metabolite.

    :param smi_out_path: The path of the file to which unique smile strings should be written representing the final
        structures generated. If None, no file is written.

    :param max_n_substructures: The maximum number of substructures to be used for building molecules.

    :param path_substructure_db: The path to the SQLite 3 substructure database, as generated by
        :py:meth:`metaboblend.databases.SubstructureDb`.

    :param path_connectivity_db: The path to the SQLite 3 connectivity database, as generated by
        :py:meth:`metaboblend.databases.create_isomorphism_database`.

    :param prescribed_mass: A mass by which to filter results; if not provided, all possible structures will be
        generated.

    :param ppm: The allowable error of the query (in parts per million). Designed to be used in accordance with
        `prescribed_mass`.

    :param out_mode: The mode in which to write to the output file.

        * **"w"** Create a new file, overwriting any existing file.

        * **"w"** Append results to an existing file.

    :param processes: How many worker processes to utilise; if left as None, :py:meth:`os.cpu_count` is used.

    :param table_name: If specified, the table specified within the substructure database will be used to generate
        molecules; else, a temporary substructure table will be created by the function.
    """

    db = SubstructureDb(path_substructure_db, path_connectivity_db)
    tolerance = 0.001
    
    if prescribed_mass is None:  # standard build method
        exact_mass__1 = round(exact_mass)
        exact_mass__0_0001 = round(exact_mass, 4)

    else:  # prescribed substructure build method
        loss = exact_mass - prescribed_mass
        exact_mass__1 = round(loss)
        exact_mass__0_0001 = round(loss, 4)

        if ((loss / 1000000) * ppm) > 0.001:
            tolerance = round((loss / 1000000) * ppm, 4)

        max_n_substructures -= 1  # we find sets of mols that add up to the loss, not the precursor mass

    if os.name == "nt":  # multiprocessing freeze support on windows
        multiprocessing.freeze_support()

    # select groups of masses at low mass resolution
    integer_mass_values = [m for m in db.select_mass_values("1", [], table_name) if m <= exact_mass__1]
    if len(integer_mass_values) == 0:
        return set()

    integer_subsets = list(subset_sum(integer_mass_values, exact_mass__1, max_n_substructures))

    configs_iso = db.k_configs(prescribed_mass is not None)
    if smi_out_path is not None:
        smi_out = open(smi_out_path, out_mode)

    substructure_subsets = []
    for integer_subset in integer_subsets:
        if len(integer_subset) > max_n_substructures or len(integer_subset) == 0:
            continue

        # refine groups of masses to 4dp mass resolution
        exact_mass_values = db.select_mass_values("0_0001", integer_subset, table_name)
        exact_subsets = []

        # use combinations to get second group of masses instead of subset sum - subset sum is integer mass only
        for mass_combo in itertools.product(*exact_mass_values):
            if abs(sum(mass_combo) - exact_mass__0_0001) <= tolerance:
                exact_subsets.append(mass_combo)

        if len(exact_subsets) == 0:
            continue

        if prescribed_mass is not None:  # add fragments mass to to loss group
            exact_subsets = [subset + (round(exact_mass - loss, 4),) for subset in exact_subsets]

        # refines groups based on ecs and gets substructures from db (appends to substructure_subsets)
        for exact_subset in exact_subsets:  
            substructure_subsets += build_from_subsets(exact_subset, mc=mc, table_name=table_name, db=db)

    with multiprocessing.Pool(processes=processes) as pool:  # send sets of substructures for building
        smi_lists = pool.map(partial(substructure_combination_build, configs_iso=configs_iso), substructure_subsets)

    smis = set([val for sublist in smi_lists for val in sublist])

    if smi_out_path is not None:
        if len(smis) != 0:
            smi_out.writelines("\n".join(smis))

        smi_out.close()

    db.close(clean)

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


def build_from_subsets(exact_subset, mc, table_name, db):
    """
    A stage of the :py:meth:`metaboblend.build_structures.build` workflow for generating molecules to a given mass
    from substructures. At this stage, mass subsets have been identified in the substructure database. Each of these
    groups are now filtered further by identifying masses that refer to valid subsets of molecules, before they are
    built to generate new molecules.

    :param db: The substructure and connectivity database. Elemental compositions and substructures are retrieved from
        the database; this information is listed as "substructure_subset" and will be appended to the 
        substructure_subsets list provided as a parameter.

    :param exact_subset: Group of masses that sum to the correct total mass, refer to substructures in the substructure
        database.

    :param mc: List of integers detailing the molecular composition of the target metabolite, in the format
        `[C, H, N, O, P, S]`.

    :param table_name: The name of the table within the substructure database from which to extract substructures. A
        prefiltered table based on the parameters specified in :py:meth:`metaboblend.build_structures.build`. See
        :py:meth:`metaboblend.build_structures.gen_subs_table`.
    """

    substructure_subsets = []

    ec_subset = combine_ecs(exact_subset, db, table_name, "0_0001")

    if len(ec_subset) == 0:
        return []

    for ec_product in itertools.product(*ec_subset):

        if mc != list(numpy.array(ec_product).sum(axis=0)):
            continue  # check each set of elemental compositions matches the target mol

        substructure_subset = db.select_substructures(ec_product, table_name)

        if len(substructure_subset) == 0:
            continue

        substructure_subsets.append(substructure_subset)

    return substructure_subsets


def substructure_combination_build(substructure_subset, configs_iso):
    """
    Final stage for building molecules; takes a combination of substructures (substructure_combination) and builds them
        according to graphs in the substructure database. May be run in parallel.

    :param substructure_subset: Combinations of substructures for building mols.

    :param configs_iso: Possible substructure combinations extracted from the connectivity database. A tuple containing
        tuples for each substructure; these tuples specify how many bonds each substructure can make.

    :return: List of smiles representing molecules generated (and the substructures used to generate them).
    """

    smis = []

    for substructure_combination in itertools.product(*substructure_subset):
        substructure_combination = sorted(substructure_combination, key=itemgetter('atoms_available', 'valence'))

        v_a = ()
        for d in substructure_combination:
            v_a += (tuple(d["degree_atoms"].values()),)  # obtain valence configuration of the set of substructures

        if str(v_a) not in configs_iso:  # check mols "fit" together according to the connectivity database
            continue

        mol_comb, atoms_available, atoms_to_remove, bond_types, bond_mismatch = reindex_atoms(substructure_combination)

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
