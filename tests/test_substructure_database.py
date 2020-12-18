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
from metaboblend.databases import *


class SubstructureDbTestCase(unittest.TestCase):
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

    def test_init(self):
        db = SubstructureDb(self.to_test_data("substructures.sqlite"),
                            self.to_test_data("connectivity.sqlite"))

        db.cursor.execute("SELECT * FROM substructures")
        first_row = db.cursor.fetchone()[0:18]
        self.assertEqual(first_row, (1, '*:c(:*)CCN', 4, 10, 56, 56.05, 56.05002399999998, 3, 6, 1, 0, 0, 0, 2,
                                     '{3: 2}', 1, '{3: [1.5, 1.5]}', '[4, 5]'))

        self.assertTrue(Chem.MolFromSmiles(first_row[1], False))
        self.assertEqual(len(db.cursor.fetchall()), 1235)

        db.cursor.execute("SELECT * FROM hmdbid_substructures")
        first_row = db.cursor.fetchone()
        self.assertEqual(first_row, ('HMDB0000073', 1))
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
        db = SubstructureDb(self.to_test_data("substructures.sqlite"))
        for i, cpd_entry in enumerate(db.select_compounds(["HMDB0000158", "HMDB0000122"])):
            self.assertLessEqual(i, 2)
            self.assertTrue(cpd_entry[0] == "HMDB0000158" or cpd_entry[0] == "HMDB0000122")

        db.close()

    def test_filter_hmdbid_substructures(self):
        db = SubstructureDb(self.to_test_data("substructures.sqlite"))
        db.filter_hmdbid_substructures(2)

        db.cursor.execute("SELECT COUNT(*) FROM filtered_hmdbid_substructures GROUP BY hmdbid")
        for i, hmdbid_count in enumerate(db.cursor.fetchall()):
            self.assertGreater(hmdbid_count[0], 1)

        self.assertEqual(i, 3)

        db.close()

    def test_generate_substructure_network(self):  # also tests get_substructure_network, get_single_edge and close
        db = SubstructureDb(self.to_test_data("substructures.sqlite"))

        self.assertEqual(db.get_single_edge([3, 4, 2]), {3: {3: None, 4: 2}, 2: {3: 1, 4: 1, 2: None}, 4: {4: None}})

        std = db.generate_substructure_network(min_node_weight=2, return_networkx=True)

        db.cursor.execute("SELECT * FROM filtered_hmdbid_substructures")
        for hmdb in db.cursor.fetchall():

            self.assertTrue(hmdb[1] in std.nodes)

        db.cursor.execute("SELECT DISTINCT substructure_id FROM filtered_hmdbid_substructures")
        self.assertEqual(len(db.cursor.fetchall()), 57)
        self.assertEqual(std.number_of_nodes(), 57)

        self.assertEqual(std.number_of_edges(), 1024)

        edge_count = []
        db.cursor.execute("SELECT * FROM substructure_graph")
        for edge in db.cursor.fetchall():
            edge_count.append(std.get_edge_data(edge[0], edge[1])["weight"])

        self.assertEqual(sum(edge_count), 2048)

        db.cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        self.assertEqual(len(db.cursor.fetchall()), 6)

        db.cursor.execute("CREATE TABLE subset_substructures AS SELECT * FROM COMPOUNDS")
        db.cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        self.assertEqual(len(db.cursor.fetchall()), 7)

        db.close()

        self.assertRaises(sqlite3.ProgrammingError, lambda: db.cursor.execute("SELECT name FROM sqlite_master WHERE type='table'"))

        db = SubstructureDb(self.to_test_data("substructures.sqlite"))
        db.cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        self.assertEqual(len(db.cursor.fetchall()), 7)

        db.close()

    def test_select_mass_values(self):
        db = SubstructureDb(self.to_test_data("substructures.sqlite"))
        ests = db.select_mass_values("1", [], None)
        exacts = db.select_mass_values("0_0001", [], None)

        self.assertEqual(len(ests), 71)
        self.assertEqual(len(exacts), 117)

        for exact in exacts:
            self.assertTrue(round(exact) in ests)

        self.assertEqual(db.select_mass_values("0_0001", [50, 64, 73], None),
                         [[50.0156], [64.0313], [73.029]])
        self.assertEqual(db.select_mass_values("0_0001", [120, 87, 87], None),
                         [[120.0423], [87.0082, 87.0446], [87.0082, 87.0446]])
        self.assertEqual(db.select_mass_values("0_0001", [50, 64, 73], None),
                         [[50.0156], [64.0313], [73.029]])
        self.assertEqual(db.select_mass_values("0_0001", [55, 80, 107], None),
                         [[55.0184, 55.0422], [80.0262, 80.05], [107.0497, 107.0735]])
        self.assertEqual(db.select_mass_values("0_0001", [63, 63, 63], None),
                         [[63.0235], [63.0235], [63.0235]])

        self.assertRaises(sqlite3.OperationalError,
                          lambda: db.select_mass_values("0_0001", [63, 63, 63], ""))
        db.close()

    def test_select_mfs(self):
        db = SubstructureDb(self.to_test_data("substructures.sqlite"))
        self.assertEqual(db.select_mfs(107.0735, None, "0_0001"), [(7, 9, 1, 0, 0, 0)])
        self.assertEqual(db.select_mfs(107.0735, None, "1"), [])
        self.assertEqual(db.select_mfs(107, None, "0_0001"), [])
        self.assertEqual(db.select_mfs(107.0735, None, "1"), [])
        self.assertEqual(db.select_mfs(107, None, "1"),
                         [(7, 7, 0, 1, 0, 0), (7, 9, 1, 0, 0, 0)])

        self.assertRaises(sqlite3.OperationalError, lambda: db.select_mfs(107.0735, "", "0_0001"))

        db.close()

    def test_k_configs(self):
        db = SubstructureDb(self.to_test_data("substructures.sqlite"),
                            self.to_test_data("connectivity.sqlite"))

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
        db = SubstructureDb(self.to_test_data("substructures.sqlite"))
        self.assertEqual(db.select_substructures([[2, 5, 0, 0, 0, 0]], None), [])
        self.assertEqual(len(db.select_substructures([[4, 4, 0, 0, 0, 0]], None)[0]), 7)
        self.assertEqual(list(db.select_substructures([[4, 4, 0, 0, 0, 0]], None)[0][0].keys()),
                         ['smiles', 'mol', 'bond_types', 'degree_atoms', 'valence', 'atoms_available', 'dummies'])

        substructures = list(db.select_substructures([[4, 4, 0, 0, 0, 0]], None)[0][0].values())
        self.assertEqual([item for i, item in enumerate(substructures) if i != 1],
                         ['*Cc(:*)cc:*',
                          {1: [1.0], 2: [1.5], 5: [1.5]},
                          {1: 1, 2: 1, 5: 1},
                          3,
                          3,
                          [0, 3, 4]])

        self.assertEqual(len(db.select_substructures([[7, 7, 0, 0, 0, 0]], None)[0]), 3)
        self.assertEqual(list(db.select_substructures([[7, 7, 0, 0, 0, 0]], None)[0][0].keys()),
                         ['smiles', 'mol', 'bond_types', 'degree_atoms', 'valence', 'atoms_available', 'dummies'])
        substructures = list(db.select_substructures([[7, 7, 0, 0, 0, 0]], None)[0][0].values())
        self.assertEqual([item for i, item in enumerate(substructures) if i != 1],
                         ['*CCc1c:*:c(*)cc1',
                          {1: [1.0], 4: [1.5], 6: [1.5, 1.0]},
                          {1: 1, 4: 1, 6: 2},
                          4,
                          3,
                          [0, 5, 7]])

        self.assertRaises(sqlite3.OperationalError,
                          lambda: db.select_substructures([[2, 5, 0, 0, 0, 0]], ""))
        db.close()

    def test_create_compound_database(self):  # also tests create_indexes
        db = SubstructureDb(self.to_test_results("substructures_new.sqlite"))
        db.create_compound_database()
        db.cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        self.assertEqual(len(db.cursor.fetchall()), 4)

        db.create_indexes()
        db.close()

        shutil.copyfile(self.to_test_data("substructures.sqlite"), self.to_test_results("substructures_copy.sqlite"))
        db = SubstructureDb(self.to_test_results("substructures_copy.sqlite"),
                            self.to_test_data("connectivity.sqlite"))
        db.create_indexes()
        db.create_compound_database()
        db.cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        self.assertEqual(len(db.cursor.fetchall()), 4)

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


if __name__ == '__main__':
    unittest.main()
