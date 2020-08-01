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
import zipfile
from metaboblend.build_structures import *
from metaboblend.databases import *


class BuildStructuresTestCase(unittest.TestCase):
    temp_results_dir = None
    temp_results_name = None

    @classmethod
    def to_test_result(cls, *args):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), cls.temp_results_name, *args)

    @classmethod
    def setUpClass(cls):
        cls.temp_results_dir = tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.realpath(__file__)))
        cls.temp_results_name = cls.temp_results_dir.name

        for compr_data in ["connectivity.zip", "test_mols.zip", "substructures.zip"]:
            zip_ref = zipfile.ZipFile(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                                   "data",
                                                   compr_data
                                                   ), 'r')
            zip_ref.extractall(cls.to_test_result())
            zip_ref.close()

    def test_build(self):
        db = SubstructureDb(self.to_test_result("substructures.sqlite"), self.to_test_result("connectivity", "pkls"),
                            self.to_test_result("connectivity", "k_graphs.sqlite"))

        smis = [{'NCCc1cc(O)ccc1O', 'NCCc1cccc(O)c1O', 'NCCc1cc(O)cc(O)c1', 'NCCc1ccc(O)c(O)c1'},
                None,
                {'N[C@@H](Cc1ccc(O)cc1)C(=O)O', 'N[C@@H](Cc1cccc(O)c1)C(=O)O', 'N[C@H](Cc1ccc(O)cc1)C(=O)O'},
                None]

        std_lens = [4, 51, 3, 1892]

        fragments = [56.05, 60.0211, 68.0262, 56.0262]

        exp_lens = [1, 41, 2, 0]

        with open(self.to_test_result("test_mols", "test_hmdbs.dictionary"), "rb") as test_hmdbs:
            record_dicts = pickle.load(test_hmdbs)
            for i, record_dict in enumerate(record_dicts.values()):
                record_dict["mol"] = Chem.MolFromSmiles(record_dict["smiles"])

                build(mc=[record_dict["C"], record_dict["H"], record_dict["N"], record_dict["O"], record_dict["P"],
                          record_dict["S"]], exact_mass=record_dict["exact_mass"],
                      smi_out_path=self.to_test_result(record_dict["HMDB_ID"] + ".smi"), max_n_substructures=3,
                      path_connectivity_db=self.to_test_result("connectivity", "k_graphs.sqlite"),
                      path_substructure_db=self.to_test_result("substructures.sqlite"),
                      prescribed_mass=None, ppm=None, out_mode="w", processes=None, table_name="substructures")

                j = 0
                unique_smis = set()
                with open(self.to_test_result(record_dict["HMDB_ID"] + ".smi"), "r") as smi_out:
                    for line in smi_out:
                        j += 1
                        unique_smis.add(line.split()[0])

                self.assertEqual(j, std_lens[i])

                if smis[i] is not None:
                    self.assertEqual(unique_smis, smis[i])
                else:
                    self.assertTrue(len(unique_smis) == 51 or len(unique_smis) == 1892)

                if os.path.isfile(self.to_test_result(record_dict["HMDB_ID"] + ".smi")):
                    os.remove(self.to_test_result(record_dict["HMDB_ID"] + ".smi"))

                build(mc=[record_dict["C"], record_dict["H"], record_dict["N"], record_dict["O"], record_dict["P"],
                          record_dict["S"]], exact_mass=record_dict["exact_mass"],
                      smi_out_path=self.to_test_result(record_dict["HMDB_ID"] + ".smi"), max_n_substructures=3,
                      prescribed_mass=fragments[i], ppm=15,
                      path_connectivity_db=self.to_test_result("connectivity", "k_graphs.sqlite"),
                      path_substructure_db=self.to_test_result("substructures.sqlite"),
                      out_mode="w", processes=None, table_name="substructures")

                j = 0
                unique_smis = set()
                with open(self.to_test_result(record_dict["HMDB_ID"] + ".smi"), "r") as smi_out:
                    for line in smi_out:
                        j += 1
                        unique_smis.add(line.split()[0])

                if i == 2:
                    self.assertEqual(unique_smis, {'N[C@@H](Cc1ccc(O)cc1)C(=O)O', 'N[C@@H](Cc1cccc(O)c1)C(=O)O'})

                self.assertEqual(len(unique_smis), exp_lens[i])

                if os.path.isfile(self.to_test_result(record_dict["HMDB_ID"] + ".smi")):
                    os.remove(self.to_test_result(record_dict["HMDB_ID"] + ".smi"))

        db.close()

    def test_gen_subs_table(self):
        db = SubstructureDb(self.to_test_result("substructures.sqlite"), "")
        table_name = gen_subs_table(db, range(5, 7), 4, 2, 500)

        i = 0
        db.cursor.execute("SELECT heavy_atoms, valence, atoms_available FROM %s" % table_name)
        for row in db.cursor.fetchall():
            i += 1
            self.assertTrue(row[0] in range(5, 7))
            self.assertTrue(row[1] <= 4)
            self.assertTrue(row[2] <= 2)

        self.assertEqual(i, 58)

        db.close()

    def test_subset_sum(self):  # also tests find_path
        self.assertEqual([s_sum for s_sum in subset_sum([1, 2, 3, 4], 5)], [[2, 3], [1, 4]])

        self.assertEqual(len(list(subset_sum(list(range(60)), 70, 3))), 378)
        self.assertEqual(len(list(subset_sum(list(range(60)), 70, 1000))), 29884)

    def test_combine_ecs(self):
        db = SubstructureDb(self.to_test_result("substructures.sqlite"), "")
        self.assertEqual(combine_ecs([54.0106, 69.0578], db, "substructures", "0_0001"),
                         [[(3, 2, 0, 1, 0, 0)], [(4, 7, 1, 0, 0, 0)]])
        self.assertEqual(combine_ecs([54, 69], db, "substructures", "1"),
                         [[(3, 2, 0, 1, 0, 0)], [(4, 7, 1, 0, 0, 0), (4, 5, 0, 1, 0, 0)]])
        self.assertEqual(combine_ecs([54.0101, 69.0580], db, "substructures", "0_0001"), [])

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
            mol_e = add_bonds(mol_comb[i], edges[i], atoms_available[i], bond_types[i])

            if i == 0:
                self.assertTrue(mol_e is None)
            else:
                self.assertEqual(Chem.MolToSmiles(mol_e.GetMol(), False), mol_out[i])


if __name__ == '__main__':
    unittest.main()
