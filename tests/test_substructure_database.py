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
import unittest
import tempfile
import shutil
import zipfile
from metaboblend.databases import *


class DatabasesTestCase(unittest.TestCase):
    temp_results_dir = None
    temp_results_name = None

    @classmethod
    def to_test_result(cls, *args):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), cls.temp_results_name, *args)

    @classmethod
    def setUpClass(cls):
        cls.temp_results_dir = tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.realpath(__file__)))
        cls.temp_results_name = cls.temp_results_dir.name

        zip_ref = zipfile.ZipFile(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                               "data",
                                               "connectivity.zip"
                                               ), 'r')
        zip_ref.extractall(cls.to_test_result())
        zip_ref.close()

        zip_ref = zipfile.ZipFile(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                               "data",
                                               "test_mols.zip"
                                               ), 'r')
        zip_ref.extractall(cls.to_test_result())
        zip_ref.close()

        zip_ref = zipfile.ZipFile(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                               "data",
                                               "substructures.zip"
                                               ), 'r')
        zip_ref.extractall(cls.to_test_result())
        zip_ref.close()

    def test_init(self):
        db = SubstructureDb(self.to_test_result("substructures.sqlite"),
                            self.to_test_result("connectivity", "k_graphs.sqlite"))

        db.cursor.execute("SELECT * FROM substructures")
        first_row = db.cursor.fetchone()[0:17]
        self.assertEqual(first_row, ('*:c(:*)CCN', 4, 10, 56, 56.05, 56.05002399999998, 3, 6, 1, 0, 0, 0, 2, '{3: 2}',
                                     1, '{3: [1.5, 1.5]}', '[4, 5]'))

        self.assertTrue(Chem.MolFromSmiles(first_row[0], False))
        self.assertEqual(len(db.cursor.fetchall()), 1235)

        db.cursor.execute("SELECT * FROM hmdbid_substructures")
        first_row = db.cursor.fetchone()
        self.assertEqual(first_row, ('HMDB0000073', '*:c(:*)CCN'))
        self.assertEqual(len(db.cursor.fetchall()), 1292)

        db.cursor.execute("SELECT * FROM compounds")
        first_row = db.cursor.fetchone()
        self.assertEqual(first_row, ('HMDB0000073', 153.078979, 'C8H11NO2', 8, 11, 1, 2, 0, 0, 'NCCC1=CC(O)=C(O)C=C1'))
        self.assertEqual(len(db.cursor.fetchall()), 3)

        db.cursor.execute("SELECT * FROM graphs.subgraphs")
        first_row = db.cursor.fetchone()
        self.assertEqual(first_row[0:9], (1, 1, b'A_', 2, '(1, 1)', '(1, 1)', '((1,), (1,))', 2, 1))
        self.assertEqual(len(db.cursor.fetchall()), 107)

        db.close()

    def test_select_compounds(self):
        db = SubstructureDb(self.to_test_result("substructures.sqlite"), "")
        for i, cpd_entry in enumerate(db.select_compounds(["HMDB0000158", "HMDB0000122"])):
            self.assertLessEqual(i, 2)
            self.assertTrue(cpd_entry[0] == "HMDB0000158" or cpd_entry[0] == "HMDB0000122")

        db.close()

    def test_filter_hmdbid_substructures(self):
        db = SubstructureDb(self.to_test_result("substructures.sqlite"), "")
        db.filter_hmdbid_substructures(2)
        db.cursor.execute("SELECT * FROM unique_hmdbid")
        self.assertEqual(db.cursor.fetchall(), [('HMDB0000073',), ('HMDB0000122',), ('HMDB0000158',), ('HMDB0000186',)])

        db.cursor.execute("SELECT * FROM filtered_hmdbid_substructures")
        for i, unique_substructure in enumerate(db.cursor.fetchall()):
            self.assertLessEqual(i, 56)
            self.assertGreaterEqual(unique_substructure[1], 2)

        db.close()

    def test_generate_substructure_network(self):  # also tests close
        db = SubstructureDb(self.to_test_result("substructures.sqlite"), "")
        std = db.generate_substructure_network(method="default", min_node_weight=2, remove_isolated=False)
        extended = db.generate_substructure_network(method="extended", min_node_weight=2, remove_isolated=False)
        parent = db.generate_substructure_network(method="parent_structure_linkage", min_node_weight=2,
                                                  remove_isolated=False)

        for s in std.nodes:
            self.assertTrue(s in extended.nodes and s in parent.nodes)

        db.cursor.execute("select * from unique_hmdbid")
        edge_count = []
        for hmdb in db.cursor.fetchall():
            self.assertTrue(hmdb[0] in parent.nodes)
            edge_count.append(len(parent.edges(hmdb[0])))

        db.cursor.execute("select distinct smiles from filtered_hmdbid_substructures")
        self.assertEqual(len(db.cursor.fetchall()), 57)
        self.assertEqual(std.number_of_nodes(), 57)
        self.assertEqual(extended.number_of_nodes(), 57)
        self.assertEqual(parent.number_of_nodes() - 4, 57)

        self.assertEqual(std.number_of_edges(), 1024)
        self.assertEqual(extended.number_of_edges(), 1024)

        self.assertEqual(parent.number_of_edges(), 114)
        self.assertEqual(sum(edge_count), 114)

        db.cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        self.assertEqual(len(db.cursor.fetchall()), 5)

        db.cursor.execute("CREATE TABLE subset_substructures AS SELECT * FROM COMPOUNDS")
        db.cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        self.assertEqual(len(db.cursor.fetchall()), 6)

        db.close()

        self.assertRaises(sqlite3.ProgrammingError, lambda: db.cursor.execute("SELECT name FROM sqlite_master WHERE type='table'"))

        db = SubstructureDb(self.to_test_result("substructures.sqlite"), "")
        db.cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        self.assertEqual(len(db.cursor.fetchall()), 3)

        db.close()

    def test_select_mass_values(self):
        db = SubstructureDb(self.to_test_result("substructures.sqlite"), "")
        ests = db.select_mass_values("1", [], "substructures")
        exacts = db.select_mass_values("0_0001", [], "substructures")

        self.assertEqual(len(ests), 71)
        self.assertEqual(len(exacts), 117)

        for exact in exacts:
            self.assertTrue(round(exact) in ests)

        self.assertEqual(db.select_mass_values("0_0001", [50, 64, 73], "substructures"),
                         [[50.0156], [64.0313], [73.029]])
        self.assertEqual(db.select_mass_values("0_0001", [120, 87, 87], "substructures"),
                         [[120.0423], [87.0082, 87.0446], [87.0082, 87.0446]])
        self.assertEqual(db.select_mass_values("0_0001", [50, 64, 73], "substructures"),
                         [[50.0156], [64.0313], [73.029]])
        self.assertEqual(db.select_mass_values("0_0001", [55, 80, 107], "substructures"),
                         [[55.0184, 55.0422], [80.0262, 80.05], [107.0497, 107.0735]])
        self.assertEqual(db.select_mass_values("0_0001", [63, 63, 63], "substructures"),
                         [[63.0235], [63.0235], [63.0235]])

        self.assertRaises(sqlite3.OperationalError,
                          lambda: db.select_mass_values("0_0001", [63, 63, 63], "substrusctures"))
        db.close()

    def test_select_ecs(self):
        db = SubstructureDb(self.to_test_result("substructures.sqlite"), "")
        self.assertEqual(db.select_ecs(107.0735, "substructures", "0_0001"), [(7, 9, 1, 0, 0, 0)])
        self.assertEqual(db.select_ecs(107.0735, "substructures", "1"), [])
        self.assertEqual(db.select_ecs(107, "substructures", "0_0001"), [])
        self.assertEqual(db.select_ecs(107.0735, "substructures", "1"), [])
        self.assertEqual(db.select_ecs(107, "substructures", "1"),
                         [(7, 9, 1, 0, 0, 0), (7, 7, 0, 1, 0, 0)])

        self.assertRaises(sqlite3.OperationalError,
                          lambda: db.select_ecs(107.0735, "substrusctures", "0_0001"))

        db.close()

    def test_k_configs(self):
        db = SubstructureDb(self.to_test_result("substructures.sqlite"),
                            self.to_test_result("connectivity", "k_graphs.sqlite"))

        k_configs = db.k_configs()
        self.assertEqual(len(k_configs), 67)
        self.assertEqual(k_configs['((1,), (1,))'], [((0, 1),)])
        self.assertEqual(k_configs['((2, 2), (2, 2), (2, 2))'],
                         [((0, 2), (0, 4), (1, 3), (1, 5), (2, 4), (3, 5)),
                          ((0, 2), (0, 5), (1, 3), (1, 4), (2, 5), (3, 4)),
                          ((0, 3), (0, 5), (1, 2), (1, 4), (2, 4), (3, 5)),
                          ((0, 3), (0, 4), (1, 2), (1, 5), (2, 5), (3, 4))])

        db.close()

    def test_select_substructures(self):
        db = SubstructureDb(self.to_test_result("substructures.sqlite"), "")
        self.assertEqual(db.select_substructures([[2, 5, 0, 0, 0, 0]], "substructures"), [])
        self.assertEqual(len(db.select_substructures([[4, 4, 0, 0, 0, 0]], "substructures")[0]), 7)
        self.assertEqual(list(db.select_substructures([[4, 4, 0, 0, 0, 0]], "substructures")[0][0].keys()),
                         ['smiles', 'mol', 'bond_types', 'degree_atoms', 'valence', 'atoms_available', 'dummies'])

        substructures = list(db.select_substructures([[4, 4, 0, 0, 0, 0]], "substructures")[0][0].values())
        self.assertEqual([item for i, item in enumerate(substructures) if i != 1],
                         ['*C(*)Cc(:*)c:*', {1: [1.0, 1.0], 3: [1.5], 6: [1.5]}, {1: 2, 3: 1, 6: 1}, 4, 3, [0, 4, 5, 7]])

        self.assertEqual(len(db.select_substructures([[7, 7, 0, 0, 0, 0]], "substructures")[0]), 3)
        self.assertEqual(list(db.select_substructures([[7, 7, 0, 0, 0, 0]], "substructures")[0][0].keys()),
                         ['smiles', 'mol', 'bond_types', 'degree_atoms', 'valence', 'atoms_available', 'dummies'])
        substructures = list(db.select_substructures([[7, 7, 0, 0, 0, 0]], "substructures")[0][0].values())
        self.assertEqual([item for i, item in enumerate(substructures) if i != 1],
                         ['*C(*)Cc1cc:*:cc1', {1: [1.0, 1.0], 5: [1.5], 7: [1.5]}, {1: 2, 5: 1, 7: 1}, 4, 3, [0, 6, 9]])

        self.assertRaises(sqlite3.OperationalError,
                          lambda: db.select_substructures([[2, 5, 0, 0, 0, 0]], "substrusctures"))
        db.close()

    def test_create_compound_database(self):  # also tests create_indexes
        db = SubstructureDb(self.to_test_result("substructures_new.sqlite"), "")
        db.create_compound_database()
        db.cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        self.assertEqual(len(db.cursor.fetchall()), 3)

        db.create_indexes()
        db.close()

        shutil.copyfile(self.to_test_result("substructures.sqlite"), self.to_test_result("substructures_copy.sqlite"))
        db = SubstructureDb(self.to_test_result("substructures_copy.sqlite"),
                            self.to_test_result("connectivity", "k_graphs.sqlite"))
        db.create_indexes()
        db.create_compound_database()
        db.cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        self.assertEqual(len(db.cursor.fetchall()), 3)

        db.cursor.execute("SELECT * FROM substructures")
        self.assertEqual(len(db.cursor.fetchall()), 0)

        db.cursor.execute("SELECT * FROM hmdbid_substructures")
        self.assertEqual(len(db.cursor.fetchall()), 0)

        db.cursor.execute("SELECT * FROM compounds")
        self.assertEqual(len(db.cursor.fetchall()), 0)

        db.cursor.execute("SELECT * FROM graphs.subgraphs")
        first_row = db.cursor.fetchone()
        self.assertEqual(first_row[0:9], (1, 1, b'A_', 2, '(1, 1)', '(1, 1)', '((1,), (1,))', 2, 1))
        self.assertEqual(len(db.cursor.fetchall()), 107)

        db.create_indexes()
        db.close()

    @classmethod
    def tearDownClass(cls):
        if cls.temp_results_dir is not None:
            cls.temp_results_dir.cleanup()


if __name__ == '__main__':
    unittest.main()
