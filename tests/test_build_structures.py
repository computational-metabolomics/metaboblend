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


import unittest
import shutil
import tempfile
from metaboblend.build_structures import *
from metaboblend.databases import *


class BuildStructuresTestCase(unittest.TestCase):
    temp_results_dir = None

    @classmethod
    def to_test_results(cls, *args):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), cls.temp_results_dir.name, *args)

    @classmethod
    def to_test_data(cls, *args):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), cls.temp_results_dir.name, "test_data", *args)

    @classmethod
    def setUpClass(cls):
        cls.temp_results_dir = tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.realpath(__file__)))

        shutil.copytree(os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_data"),
                        cls.to_test_results("test_data"))

    def test_build(self):  # core - all other build functions rely on
        db = SubstructureDb(self.to_test_data("substructures.sqlite"), self.to_test_data("connectivity.sqlite"))

        # ref data
        smis = [{'NCCc1cc(O)ccc1O', 'NCCc1cccc(O)c1O', 'NCCc1cc(O)cc(O)c1', 'NCCc1ccc(O)c(O)c1'},
                None,
                {'N[C@@H](Cc1ccc(O)cc1)C(=O)O', 'N[C@@H](Cc1cccc(O)c1)C(=O)O', 'N[C@H](Cc1ccc(O)cc1)C(=O)O'},
                None]
        std_lens = [4, 47, 3, 1892]
        fragments = [56.05, 60.0211, 68.0262, 56.0262]
        exp_lens = [1, 41, 2, 0]

        # hmdb records to build from
        with open(self.to_test_data("test_hmdbs.dictionary"), "rb") as test_hmdbs:
            record_dicts = pickle.load(test_hmdbs)

            for i, record_dict in enumerate(record_dicts.values()):

                # test standard building
                built_smis = build(
                    mf=[record_dict["C"], record_dict["H"], record_dict["N"],
                        record_dict["O"], record_dict["P"], record_dict["S"]],
                    exact_mass=record_dict["exact_mass"], max_n_substructures=3, db=db, ppm=None, ncpus=None,
                    table_name=None, isomeric_smiles=True, retain_substructures=True, prescribed_substructures=None,
                    tolerance=None
                )

                self.assertEqual(len(built_smis), std_lens[i])

                if smis[i] is not None:
                    self.assertEqual(set(built_smis.keys()), smis[i])

                # test prescribed substructure building
                built_smis = build(
                    mf=[record_dict["C"], record_dict["H"], record_dict["N"],
                        record_dict["O"], record_dict["P"], record_dict["S"]],
                    exact_mass=record_dict["exact_mass"], max_n_substructures=3, tolerance=0.0001,
                    prescribed_substructures=get_possible_fragment_ions(fragments[i], db, 2, 5, 0.0001), ppm=15,
                    ncpus=None, isomeric_smiles=True, retain_substructures=False, db=db, table_name=None
                )

                if i == 2:
                    self.assertEqual(set(built_smis.keys()), {'N[C@@H](Cc1ccc(O)cc1)C(=O)O', 'N[C@@H](Cc1cccc(O)c1)C(=O)O'})

                self.assertEqual(len(built_smis.keys()), exp_lens[i])

        db.close()

    def test_substructure_combination_build(self):
        db = SubstructureDb(self.to_test_data("substructures.sqlite"),
                            self.to_test_data("connectivity.sqlite"))

        ec_products = [[(3, 3, 0, 1, 0, 0), (5, 8, 1, 1, 0, 0)],
                       [(3, 5, 0, 1, 0, 0), (2, 4, 1, 1, 0, 0), (4, 2, 0, 1, 0, 0)],
                       [(3, 5, 1, 0, 0, 0), (2, 3, 0, 2, 0, 0)],
                       [(3, 6, 0, 4, 0, 0), (5, 9, 0, 3, 0, 0), (4, 7, 0, 4, 0, 0)]]
        configs_iso = db.k_configs()
        lens = [3, 0, 0, 44]

        for i, ec_product in enumerate(ec_products):
            substructure_subset = db.select_substructures(ec_product, None)

            smis = substructure_combination_build(substructure_subset, configs_iso, prescribed_method=False,
                                                  isomeric_smiles=True, bond_enthalpies=get_bond_enthalpies(),
                                                  retain_substructures=False)

            self.assertEqual(len(smis.keys()), lens[i])

            if i == 0:
                self.assertEqual(list(smis.keys()), ['NCCc1ccc(O)c(O)c1', 'NCCc1cc(O)ccc1O', 'NCCc1cc(O)cc(O)c1'])

        db.close()

    def test_build_from_subsets(self):
        db = SubstructureDb(self.to_test_data("substructures.sqlite"))

        mcs = [[12, 22, 0, 11, 0, 0], [10, 0, 0, 0, 0, 0],
               [9, 11, 1, 3, 0, 0], [9, 11, 1, 3, 0, 0]]
        exact_subsets = [(103.0395, 119.0344, 120.0423),
                         (84.0449, 97.029), (50.0156, 57.0215, 74.0368), (50.0156, 57.034, 74.0242)]

        lens = [13, 0, 1, 1]
        
        for i, mc, exact_subset in zip(range(len(mcs)), mcs, exact_subsets):
            substructure_subsets = build_from_subsets(exact_subset, mc, None, db)

            if i == 1:
                self.assertEqual(len(substructure_subsets), 0)
            else:
                self.assertEqual(len(substructure_subsets[0][0]), lens[i])

            if i == 2:
                del substructure_subsets[0][0][0]["mol"]
                self.assertEqual(substructure_subsets[0][0],
                                 [{'atoms_available': 2,
                                   'bond_types': {1: [1.0, 1.5], 4: [1.5, 1.0]},
                                   'degree_atoms': {1: 2, 4: 2},
                                   'dummies': [0, 2, 3, 5],
                                   'smiles': '*c1:*:*:c(*)cc1',
                                   'valence': 4}])

        db.close()

    def test_generate_structures(self):  # tests vs build
        db = SubstructureDb(self.to_test_data("substructures.sqlite"), self.to_test_data("connectivity.sqlite"))

        fragments = [56.05, 60.0211, 68.0262, 56.0262]

        with open(self.to_test_data("test_hmdbs.dictionary"), "rb") as test_hmdbs:
            record_dicts = pickle.load(test_hmdbs)
            for i, record_dict in enumerate(record_dicts.values()):
                ms_data = {record_dict["HMDB_ID"]: {"mf": [record_dict["C"], record_dict["H"], record_dict["N"],
                                                           record_dict["O"], record_dict["P"], record_dict["S"]],
                                                    "exact_mass": record_dict["exact_mass"]}}

                # test standard building
                returned_smis = list(
                    generate_structures(ms_data, path_substructure_db=self.to_test_data("substructures.sqlite"),
                                        write_csv_output=True, path_out=self.to_test_results(),
                                        max_degree=6, max_atoms_available=2, max_n_substructures=3,
                                        path_connectivity_db=self.to_test_data("connectivity.sqlite"),
                                        minimum_frequency=None, yield_smis=True, isomeric_smiles=True,
                                        retain_substructures=True))

                returned_smis = returned_smis[0][record_dict["HMDB_ID"]]

                build_smis = build(
                    mf=[record_dict["C"], record_dict["H"], record_dict["N"],
                        record_dict["O"], record_dict["P"], record_dict["S"]],
                    exact_mass=record_dict["exact_mass"],
                    max_n_substructures=3, ppm=None, ncpus=None, table_name=None, isomeric_smiles=True,
                    retain_substructures=True, db=db, tolerance=0.0001, prescribed_substructures=None
                )

                self.assertEqual(set(build_smis.keys()), set(returned_smis))

                ms_data = {record_dict["HMDB_ID"]: {"mf": [record_dict["C"], record_dict["H"], record_dict["N"],
                                                           record_dict["O"], record_dict["P"], record_dict["S"]],
                                                    "exact_mass": record_dict["exact_mass"],
                                                    "prescribed_mass": fragments[i]}}

                # test prescribed building
                returned_smis = list(
                    generate_structures(ms_data, path_substructure_db=self.to_test_data("substructures.sqlite"),
                                        write_csv_output=True, path_out=self.to_test_results(),
                                        max_degree=6, max_atoms_available=2, max_n_substructures=3,
                                        path_connectivity_db=self.to_test_data("connectivity.sqlite"),
                                        minimum_frequency=None, yield_smis=True, isomeric_smiles=True,
                                        retain_substructures=False))

                returned_smis = returned_smis[0][record_dict["HMDB_ID"]]

                prescribed_substructures = get_possible_fragment_ions(fragments[i], db)

                build_smis = build(
                    mf=[record_dict["C"], record_dict["H"], record_dict["N"],
                        record_dict["O"], record_dict["P"], record_dict["S"]],
                    exact_mass=record_dict["exact_mass"], max_n_substructures=3, ppm=0, retain_substructures=False,
                    ncpus=None, table_name=None, isomeric_smiles=True, db=db, tolerance=0.0001,
                    prescribed_substructures=prescribed_substructures
                )

                self.assertEqual(set(build_smis.keys()), set(returned_smis))

            ms_data = {}
            for i, record_dict in enumerate(record_dicts.values()):
                record_dict["mol"] = Chem.MolFromSmiles(record_dict["smiles"])
                ms_data[record_dict["HMDB_ID"]] = {"mf": [record_dict["C"], record_dict["H"], record_dict["N"],
                                                          record_dict["O"], record_dict["P"], record_dict["S"]],
                                                   "exact_mass": record_dict["exact_mass"],
                                                   "prescribed_masses": None}

            # test building with multiple inputs
            returned_smi_list = list(
                generate_structures(ms_data, path_substructure_db=self.to_test_data("substructures.sqlite"),
                                    write_csv_output=True, path_out=self.to_test_results(),
                                    max_degree=6, max_atoms_available=2, max_n_substructures=3,
                                    path_connectivity_db=self.to_test_data("connectivity.sqlite"),
                                    minimum_frequency=None, yield_smis=True, isomeric_smiles=True,
                                    retain_substructures=False))

            for i, record_dict in enumerate(record_dicts.values()):
                build_smis = build(
                    mf=[record_dict["C"], record_dict["H"], record_dict["N"],
                        record_dict["O"], record_dict["P"], record_dict["S"]],
                    exact_mass=record_dict["exact_mass"], retain_substructures=False,
                    max_n_substructures=3, ppm=None, ncpus=None, table_name=None, isomeric_smiles=True,
                    prescribed_substructures=None, db=db, tolerance=0.0001
                )

                self.assertEqual(set(build_smis.keys()), set(returned_smi_list[i][record_dict["HMDB_ID"]]))

        db.close()

    def test_annotate_msn(self):  # tests vs build_msn
        db = SubstructureDb(self.to_test_data("substructures.sqlite"))

        overall_lens = [3, 41, 2, 0]
        smis = [{'NCCc1ccc(O)c(O)c1', 'NCCc1cc(O)ccc1O', 'NCCc1cc(O)cc(O)c1'},
                None,
                {'N[C@@H](Cc1cccc(O)c1)C(=O)O', 'N[C@@H](Cc1ccc(O)cc1)C(=O)O'},
                set()]
        freqs = [1, 0, 0, 0]

        fragments = [56.05, 60.0211, 68.0262, 56.0262]

        with open(self.to_test_data("test_hmdbs.dictionary"), "rb") as test_hmdbs:
            record_dicts = pickle.load(test_hmdbs)
            for i, record_dict in enumerate(record_dicts.values()):

                if not os.path.exists(self.to_test_results("annotate")):
                    os.mkdir(self.to_test_results("annotate"))

                ms_data = {record_dict["HMDB_ID"]: {"mf": [record_dict["C"], record_dict["H"], record_dict["N"],
                                                           record_dict["O"], record_dict["P"], record_dict["S"]],
                                                    "exact_mass": record_dict["exact_mass"],
                                                    "neutral_fragment_masses": fragments}}

                # test standard building
                returned_smis = list(annotate_msn(
                    ms_data, max_degree=6, max_atoms_available=2, max_n_substructures=3,
                    write_csv_output=True, retain_substructures=False, path_out=self.to_test_results(),
                    path_connectivity_db=self.to_test_data("connectivity.sqlite"),
                    path_substructure_db=self.to_test_data("substructures.sqlite"),
                    minimum_frequency=None, yield_smis=True, isomeric_smiles=True
                ))

                returned_smis = returned_smis[0][record_dict["HMDB_ID"]]

                self.assertEqual(len([t[1] for t in returned_smis if t[1] > 1]), freqs[i])

                if smis[i] is not None:
                    self.assertEqual(set(t[0] for t in returned_smis), smis[i])

                if i == 0:
                    self.assertEqual(returned_smis[2][1], 3)

            ms_data = {}
            for i, record_dict in enumerate(record_dicts.values()):
                record_dict["mol"] = Chem.MolFromSmiles(record_dict["smiles"])
                ms_data[record_dict["HMDB_ID"]] = {"mf": [record_dict["C"], record_dict["H"], record_dict["N"],
                                                          record_dict["O"], record_dict["P"], record_dict["S"]],
                                                   "exact_mass": record_dict["exact_mass"],
                                                   "neutral_fragment_masses": fragments}

            os.mkdir(self.to_test_results("annotate_multi"))

            # test building with multiple inputs
            returned_smi_list = list(annotate_msn(
                ms_data, max_degree=6, max_atoms_available=2, max_n_substructures=3,
                path_out=self.to_test_results("annotate_multi"), write_csv_output=True,
                path_connectivity_db=self.to_test_data("connectivity.sqlite"),
                path_substructure_db=self.to_test_data("substructures.sqlite"),
                minimum_frequency=None, yield_smis=True,
                isomeric_smiles=True, retain_substructures=False
            ))

            for i, record_dict in enumerate(record_dicts.values()):
                self.assertEqual(len(returned_smi_list[i][record_dict["HMDB_ID"]]), overall_lens[i])
        
        db.close()

    def test_results_db(self):
        fragments = [56.05, 60.0211, 68.0262, 56.0262]

        with open(self.to_test_data("test_hmdbs.dictionary"), "rb") as test_hmdbs:
            record_dicts = pickle.load(test_hmdbs)

        ms_data = {}
        for i, record_dict in enumerate(record_dicts.values()):
            record_dict["mol"] = Chem.MolFromSmiles(record_dict["smiles"])
            ms_data[record_dict["HMDB_ID"]] = {"mf": [record_dict["C"], record_dict["H"], record_dict["N"],
                                                      record_dict["O"], record_dict["P"], record_dict["S"]],
                                               "exact_mass": record_dict["exact_mass"],
                                               "neutral_fragment_masses": fragments}

        os.mkdir(self.to_test_results("test_results_db"))

        list(annotate_msn(
            ms_data, max_degree=6, max_atoms_available=2, max_n_substructures=3,
            path_out=self.to_test_results("test_results_db"), write_csv_output=True,
            path_connectivity_db=self.to_test_data("connectivity.sqlite"),
            path_substructure_db=self.to_test_data("substructures.sqlite"),
            minimum_frequency=None, yield_smis=True,
            isomeric_smiles=True, retain_substructures=True
        ))

        # is the sqlite database the size we expect?
        self.assertEqual(os.path.getsize(self.to_test_results("test_results_db", "metaboblend_results.sqlite")), 53248)

        # are the csv files the same as the reference?
        with open(self.to_test_results("test_results_db", "metaboblend_queries.csv"), "r") as results_file, \
             open(self.to_test_data("metaboblend_queries.csv"), "r") as test_file:

            for results_line, test_line in zip(results_file, test_file):
                self.assertEqual(results_line, test_line)

        with open(self.to_test_results("test_results_db", "metaboblend_structures.csv"), "r") as results_file, \
             open(self.to_test_data("metaboblend_structures.csv"), "r") as test_file:

            for results_line, test_line in zip(results_file, test_file):
                self.assertEqual(results_line, test_line)

    def test_gen_subs_table(self):
        db = SubstructureDb(self.to_test_data("substructures.sqlite"), "")
        table_name = gen_subs_table(db, 5, 6, 4, 2, 500)

        i = 0
        db.cursor.execute("SELECT heavy_atoms, valence, atoms_available FROM %s" % table_name + "_substructures")
        for row in db.cursor.fetchall():
            i += 1

            self.assertTrue(row[0] in range(5, 7))
            self.assertTrue(row[1] <= 4)
            self.assertTrue(row[2] <= 2)

        self.assertEqual(i, 57)

        db.close()

    def test_subset_sum(self):  # also tests find_path
        self.assertEqual([s_sum for s_sum in subset_sum([1, 2, 3, 4], 5)], [[2, 3], [1, 4]])

        self.assertEqual(len(list(subset_sum(list(range(60)), 70, 3))), 378)
        self.assertEqual(len(list(subset_sum(list(range(60)), 70, 1000))), 29884)

    def test_combine_ecs(self):
        db = SubstructureDb(self.to_test_data("substructures.sqlite"), "")
        self.assertEqual(combine_mfs([54.0106, 69.0578], db, None, "0_0001"),
                         [[(3, 2, 0, 1, 0, 0)], [(4, 7, 1, 0, 0, 0)]])
        self.assertEqual(combine_mfs([54, 69], db, None, "1"),
                         [[(3, 2, 0, 1, 0, 0)], [(4, 7, 1, 0, 0, 0)]])
        self.assertEqual(combine_mfs([54.0101, 69.0580], db, None, "0_0001"), [])

        db.close()

    def test_reindex_atoms(self):
        substructure_combinations = [
            [{'smiles': '*C(*)C(=O)O', 'mol': None, 'bond_types': {1: [1.0, 1.0]}, 'degree_atoms': {1: 2},
              'valence': 2, 'atoms_available': 1, 'dummies': [0, 2]},
             {'smiles': 'NCCc1c:*:*:cc1', 'mol': None, 'bond_types': {4: [1.5], 6: [1.5], 7: [1.5]},
              'degree_atoms': {4: 1, 7: 1}, 'valence': 2, 'atoms_available': 2, 'dummies': [5, 6]}],
            [{'smiles': '*[C@@H](O)[C@@H](*)O', 'mol': None, 'bond_types': {1: [1.0], 3: [1.0]},
              'degree_atoms': {1: 1, 3: 1}, 'valence': 2, 'atoms_available': 2, 'dummies': [0, 5]},
             {'smiles': 'OC1**[C@@H](O)[C@H](O)[C@H]1O', 'mol': None, 'bond_types': {0: [1.0], 3: [1.0], 4: [1.0]},
              'degree_atoms': {0: 1, 4: 1}, 'valence': 2, 'atoms_available': 2, 'dummies': [2, 3]}],
            [{'smiles': '*C[C@H](N)C(=O)O', 'mol': None, 'bond_types': {2: [1.0]}, 'degree_atoms': {2: 1},
              'valence': 1, 'atoms_available': 1, 'dummies': [3]},
             {'smiles': '*c1ccc(O)cc1', 'mol': None, 'bond_types': {1: [1.0]}, 'degree_atoms': {1: 1},
              'valence': 1, 'atoms_available': 1, 'dummies': [0]}]
        ]

        reindexed = [
            ["*C(*)C(=O)O.NCCc1c:*:*:cc1", [1, 10, 13], [0, 2, 11, 12],
             {1: [1.0, 1.0], 10: [1.5], 12: [1.5], 13: [1.5]}],
            ["*[C@@H](O)[C@@H](*)O.OC1**[C@@H](O)[C@H](O)[C@H]1O", [1, 3, 6, 10], [0, 5, 8, 9],
             {1: [1.0], 3: [1.0], 6: [1.0], 9: [1.0], 10: [1.0]}],
            ["*C[C@H](N)C(=O)O.*c1ccc(O)cc1", [2, 8], [3, 7], {2: [1.0], 8: [1.0]}]
        ]

        for substructure_combination, reindex in zip(substructure_combinations, reindexed):
            substructure_combination[0]["mol"] = Chem.MolFromSmiles(substructure_combination[0]["smiles"], False)
            substructure_combination[1]["mol"] = Chem.MolFromSmiles(substructure_combination[1]["smiles"], False)

            mol_comb, atoms_available, atoms_to_remove, bond_types, bond_mismatch = reindex_atoms(substructure_combination)
            self.assertEqual([Chem.MolToSmiles(mol_comb), atoms_available, atoms_to_remove, bond_types], reindex)

    def test_add_bonds(self):
        mol_comb = [Chem.MolFromSmiles("*C(*)C(=O)O.NCCc1c:*:*:cc1", False),
                    Chem.MolFromSmiles("*[C@@H](O)[C@@H](*)O.OC1**[C@@H](O)[C@H](O)[C@H]1O", False),
                    Chem.MolFromSmiles("*C[C@H](N)C(=O)O.*c1ccc(O)cc1", False)]

        atoms_available = [[1, 10, 13], [1, 3, 6, 10], [2, 8]]

        bond_types = [{1: [1.0, 1.0], 10: [1.5], 12: [1.5], 13: [1.5]},
                       {1: [1.0], 3: [1.0], 6: [1.0], 9: [1.0], 10: [1.0]},
                       {2: [1.0], 8: [1.0]}]

        edges = [((0, 1), (0, 2)), ((0, 2), (1, 3)), ((0, 1),)]

        mol_out = [None, "*[CH]1(O)OC2**[CH](O)(C(O)C2O)[CH]1(*)O", "*C[CH](N)(C(=O)O)c1(*)ccc(O)cc1"]

        for i in range(len(atoms_available)):
            mol_e, total_bde = add_bonds(mol_comb[i], edges[i], atoms_available[i], bond_types[i], get_bond_enthalpies())

            if i == 0:
                self.assertTrue(mol_e is None)
            else:
                self.assertEqual(Chem.MolToSmiles(mol_e.GetMol(), False), mol_out[i])


if __name__ == '__main__':
    unittest.main()
