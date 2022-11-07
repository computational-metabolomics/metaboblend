#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019-2020 Jack Gisby, Ralf Weber
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

import shutil
import tempfile
import unittest

from metaboblend.parse import reformat_xml
from metaboblend.databases.substructures import *


class DatabasesTestCase(unittest.TestCase):
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

        shutil.copytree(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "test_data"),
                        cls.to_test_results("test_data"))

        cls.maxDiff = None

    def test_reformat_xml(self):

        reformat_xml(self.to_test_data("HMDB0000073_raw.xml"))

        with open(self.to_test_data("HMDB0000073_raw.xml"), "r", encoding="utf-8") as fn_hmdb:
            xml_contents = fn_hmdb.readlines()

            self.assertTrue("hmdb" in xml_contents[1])
            self.assertTrue(len(xml_contents), 3283)

    def test_parse_xml(self):

        hmdbs = ["HMDB0000073", "HMDB0000122", "HMDB0000158", "HMDB0000186"]
        lengths = [111, 109, 112, 108]
        elements = ["C8H11NO2", "C6H12O6", "C9H11NO3", "C12H22O11"]
        smis = ["NCCC1=CC(O)=C(O)C=C1",
                "[H]C1(O)O[C@]([H])(CO)[C@@]([H])(O)[C@]([H])(O)[C@@]1([H])O",
                "N[C@@H](CC1=CC=C(O)C=C1)C(O)=O",
                "OC[C@H]1O[C@@H](O[C@H]2[C@H](O)[C@@H](O)[C@@H](O)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@H]1O"]

        with open(self.to_test_data("parsed_records.dictionary"), "rb") as parsed:
            parsed_records = pickle.load(parsed)

        for i, hmdb in enumerate(["HMDB0000073", "HMDB0000122", "HMDB0000158", "HMDB0000186"]):

            for record_out in parse_xml(self.to_test_data(hmdb + ".xml")):

                self.assertEqual(len(record_out), lengths[i])
                self.assertEqual(record_out["accession"], hmdbs[i])
                self.assertEqual(record_out["smiles"], smis[i])
                self.assertEqual(record_out["chemical_formula"], elements[i])
                self.assertEqual(record_out, parsed_records[hmdb + ".xml"])

    def test_filter_records(self):

        with open(self.to_test_data("parsed_records.dictionary"), "rb") as p:
            parsed_records = pickle.load(p)

        with open(self.to_test_data("test_hmdbs.dictionary"), "rb") as test_hmdbs:
            filtered_records = pickle.load(test_hmdbs)

        record_gen = filter_records(parsed_records.values(), isomeric_smiles=True)
        test_filtered_records = {}
        for record in record_gen:

            # don't check smiles
            for field in ["mol", "smiles_rdkit", "smiles_rdkit_kek"]:
                del record[field]

            test_filtered_records[record["HMDB_ID"]] = record

        self.assertEqual(test_filtered_records, filtered_records)

    def test_get_substructure_bond_idx(self):

        with open(self.to_test_data("test_hmdbs.dictionary"), "rb") as test_hmdbs:
            record_dict = pickle.load(test_hmdbs)["HMDB0000186"]
            mol = Chem.MolFromSmiles(record_dict["smiles"])

        subs_mol = Chem.MolFromSmiles("OC[C@H]1OC[C@H](O)[C@@H](O)[C@H]1O")
        self.assertEqual(get_substructure_bond_idx(subs_mol, mol), (0, 1, 2, 22, 3, 16, 17, 18, 19, 20, 21))

        subs_mol = Chem.MolFromSmiles("OC[C@@H]1C[C@H](O)[C@@H](O)[C@@H](O)O1")
        self.assertEqual(get_substructure_bond_idx(subs_mol, mol), (0, 1, 2, 22, 3, 4, 16, 17, 18, 19, 20))

    def test_subset_sgs_sizes(self):

        sgs = [[(0, 1, 2, 22, 3, 16, 17, 18, 19, 20, 21), (0, 1, 2, 22, 3, 4, 16, 17, 18, 19, 20)]]

        self.assertEqual(len(subset_sgs_sizes(sgs, -10, 100)[0]), 2)
        self.assertEqual(len(subset_sgs_sizes(sgs, 100, 100)), 0)
        self.assertEqual(len(subset_sgs_sizes(sgs, 0, 0)), 0)
        self.assertEqual(len(subset_sgs_sizes(sgs, 11, 11)[0]), 2)
        self.assertEqual(len(subset_sgs_sizes(sgs, 11, 12)[0]), 2)
        self.assertEqual(len(subset_sgs_sizes(sgs, 10, 11)[0]), 2)
        self.assertEqual(len(subset_sgs_sizes(sgs, 12, 100)), 0)
        self.assertEqual(len(subset_sgs_sizes(sgs, 0, 10)), 0)

    def test_get_sgs(self):

        with open(self.to_test_data("test_hmdbs.dictionary"), "rb") as test_hmdbs:
            record_dict = pickle.load(test_hmdbs)["HMDB0000186"]
            record_dict["mol"] = Chem.MolFromSmiles(record_dict["smiles"])
            mol_ids = [bond.GetIdx() for bond in record_dict["mol"].GetBonds()]

        sgs = get_sgs(record_dict, 2, 9, method="exhaustive")
        for edges in sgs:
            for edge_set in edges:
                self.assertTrue(2 <= len(edge_set) <= 9)
                [self.assertTrue(bond in mol_ids) for bond in edge_set]

        sgs = get_sgs(record_dict, 0, 20, method="RECAP")
        for edges in sgs:
            for edge_set in edges:
                self.assertTrue(0 <= len(edge_set) <= 20)
                [self.assertTrue(bond in mol_ids) for bond in edge_set]

        sgs = get_sgs(record_dict, 0, 20, method="BRICS")
        for edges in sgs:
            for edge_set in edges:
                self.assertTrue(0 <= len(edge_set) <= 20)
                [self.assertTrue(bond in mol_ids) for bond in edge_set]

    def test_get_substructure(self):

        with open(self.to_test_data("test_hmdbs.dictionary"), "rb") as test_hmdbs:
            record_dict = pickle.load(test_hmdbs)["HMDB0000186"]
            mol = Chem.MolFromSmiles(record_dict["smiles"])

        libs = [{'smiles': '*[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O', 'bond_types': {4: [1.0]},
                 'degree_atoms': {4: 1}, 'valence': 1, 'atoms_available': 1, 'dummies': [5]},
                {'smiles': '*O[C@@H]1O[C@H](CO)[C@H](*)[C@H](O)[C@H]1O', 'bond_types': {5: [1.0], 11: [1.0]},
                 'degree_atoms': {5: 1, 11: 1}, 'valence': 2, 'atoms_available': 2, 'dummies': [6, 12]}]

        for edges in [[(0, 1, 2, 22, 3, 16, 17, 18, 19, 20, 21), (0, 1, 2, 22, 3, 4, 16, 17, 18, 19, 20)]]:
            for i, edge_idx in enumerate(edges):
                lib = get_substructure(mol, edge_idx, isomeric_smiles=True)
                del lib["mol"]
                self.assertEqual(lib, libs[i])

    def test_get_elements(self):

        compositions = [{'C': 8, 'H': 11, 'N': 1, 'O': 2, 'P': 0, 'S': 0, '*': 0},
                        {'C': 6, 'H': 12, 'N': 0, 'O': 6, 'P': 0, 'S': 0, '*': 0},
                        {'C': 9, 'H': 11, 'N': 1, 'O': 3, 'P': 0, 'S': 0, '*': 0},
                        {'C': 12, 'H': 22, 'N': 0, 'O': 11, 'P': 0, 'S': 0, '*': 0}]

        with open(self.to_test_data("test_hmdbs.dictionary"), "rb") as test_hmdbs:
            record_dicts = pickle.load(test_hmdbs)
            for i, record_dict in enumerate(record_dicts.values()):
                mol = Chem.MolFromSmiles(record_dict["smiles"])

                self.assertEqual(get_elements(mol), compositions[i])

    def test_calculate_exact_mass(self):

        masses = [153.07897899999998, 180.06338999999997, 181.07389399999997, 342.1162150000005]

        with open(self.to_test_data("test_hmdbs.dictionary"), "rb") as test_hmdbs:
            record_dicts = pickle.load(test_hmdbs)
            for i, record_dict in enumerate(record_dicts.values()):
                mol = Chem.MolFromSmiles(record_dict["smiles"])

                self.assertEqual(calculate_exact_mass(mol), masses[i])

        ref_db = sqlite3.connect(self.to_test_data("substructures.sqlite"))
        ref_db_cursor = ref_db.cursor()
        ref_db_cursor.execute("SELECT exact_mass__0_0001, mol FROM substructures")
        for row in ref_db_cursor.fetchall():
            self.assertEqual(round(calculate_exact_mass(Chem.Mol(row[1])), 4), row[0])

        ref_db.close()

    def test_create_substructure_database(self):

        records = [self.to_test_data(r + ".xml") for r in ["HMDB0000073", "HMDB0000122", "HMDB0000158", "HMDB0000186"]]

        create_substructure_database(records, self.to_test_results("test_db.sqlite"), 4, 8, method="exhaustive", isomeric_smiles=True)
        test_db = sqlite3.connect(self.to_test_results("test_db.sqlite"))

        test_db_cursor = test_db.cursor()

        test_db_cursor.execute("""SELECT smiles,
                                         heavy_atoms,
                                         length,
                                         exact_mass__1,
                                         exact_mass__0_0001,
                                         exact_mass,
                                         C,
                                         H,
                                         N,
                                         O,
                                         P,
                                         S,
                                         valence,
                                         valence_atoms,
                                         atoms_available,
                                         bond_types,
                                         dummies
                                  FROM substructures WHERE valence <= 4""")

        for i, row in enumerate(test_db_cursor.fetchall()):
            if i == 0:
                self.assertEqual(row, ('*:c(:*)CCN', 4, 10, 56, 56.05, 56.05002399999998, 3, 6, 1, 0, 0, 0, 2, '{3: 2}',
                                       1, '{3: [1.5, 1.5]}', '[4, 5]')
                                 )

            total_rows = i

        self.assertEqual(total_rows, 139)

        test_db_cursor.execute("SELECT * FROM hmdbid_substructures")
        for i, row in enumerate(test_db_cursor.fetchall()):
            if i == 0:
                self.assertEqual(row, ('HMDB0000073', 1))
            total_rows = i

        self.assertEqual(total_rows, 150)

        test_db_cursor.execute("SELECT * FROM compounds")
        for i, row in enumerate(test_db_cursor.fetchall()):
            if i == 0:
                self.assertEqual(row,
                                 ('HMDB0000073', 153.078979, 'C8H11NO2', 8, 11, 1, 2, 0, 0, 'NCCC1=CC(O)=C(O)C=C1'))
            total_rows = i

        self.assertEqual(total_rows, 3)

        test_db_cursor.execute("SELECT heavy_atoms FROM substructures")
        unique_ha = set()
        for ha in test_db_cursor.fetchall():
            self.assertTrue(4 <= ha[0] <= 8)
            unique_ha.add(ha[0])

        [self.assertTrue(ha in unique_ha) for ha in [4, 5, 6, 7, 8]]

        test_db.close()

    def test_update_substructure_database(self):  # requires create_compound_database from SubstructureDb

        db = SubstructureDb(self.to_test_results("test_db.sqlite"), "")
        db.create_compound_database()
        db.close()

        for record in ["HMDB0000073", "HMDB0000122", "HMDB0000158", "HMDB0000186"]:
            record = self.to_test_data(record + ".xml")

            update_substructure_database(self.to_test_data(record),
                                         self.to_test_results("test_db.sqlite"), 4, 8,
                                         method="exhaustive", isomeric_smiles=True)

        shutil.copyfile(self.to_test_data("substructures.sqlite"), self.to_test_results("substructures_copy.sqlite"))

        test_db = sqlite3.connect(self.to_test_results("test_db.sqlite"))
        test_db_cursor = test_db.cursor()

        test_db_cursor.execute("""SELECT smiles,
                                         heavy_atoms,
                                         length,
                                         exact_mass__1,
                                         exact_mass__0_0001,
                                         exact_mass,
                                         C,
                                         H,
                                         N,
                                         O,
                                         P,
                                         S,
                                         valence,
                                         valence_atoms,
                                         atoms_available,
                                         bond_types,
                                         dummies
                                  FROM substructures WHERE valence <= 4""")

        for i, row in enumerate(test_db_cursor.fetchall()):
            if i == 0:
                self.assertEqual(row, ('*:c(:*)CCN', 4, 10, 56, 56.05, 56.05002399999998, 3, 6, 1, 0, 0, 0, 2, '{3: 2}',
                                       1, '{3: [1.5, 1.5]}', '[4, 5]'))

            total_rows = i

        self.assertEqual(total_rows, 575)

        test_db_cursor.execute("SELECT * FROM hmdbid_substructures")
        for i, row in enumerate(test_db_cursor.fetchall()):
            if i == 0:
                self.assertEqual(row, ('HMDB0000073', 1))
            total_rows = i

        self.assertEqual(total_rows, 1292)

        test_db_cursor.execute("SELECT * FROM compounds")
        for i, row in enumerate(test_db_cursor.fetchall()):
            if i == 0:
                self.assertEqual(row,
                                 ('HMDB0000073', 153.078979, 'C8H11NO2', 8, 11, 1, 2, 0, 0, 'NCCC1=CC(O)=C(O)C=C1'))
            total_rows = i

        self.assertEqual(total_rows, 3)

        test_db_cursor.execute("SELECT heavy_atoms FROM substructures")
        unique_ha = set()
        for ha in test_db_cursor.fetchall():
            self.assertTrue(4 <= ha[0] <= 8)
            unique_ha.add(ha[0])

        [self.assertTrue(ha in unique_ha) for ha in [4, 5, 6, 7, 8]]

        test_db.close()

        # small substructures
        db = SubstructureDb(self.to_test_results("test_db.sqlite"), "")
        db.create_compound_database()
        db.close()

        for record in ["HMDB0000073", "HMDB0000122", "HMDB0000158", "HMDB0000186"]:
            record = self.to_test_data(record + ".xml")

            update_substructure_database(self.to_test_data(record),
                                         self.to_test_results("test_db.sqlite"), 1, 1,
                                         method="exhaustive", isomeric_smiles=True)

        test_db = sqlite3.connect(self.to_test_results("test_db.sqlite"))
        test_db_cursor = test_db.cursor()

        test_db_cursor.execute("""SELECT smiles,
                                                 heavy_atoms,
                                                 length,
                                                 exact_mass__1,
                                                 exact_mass__0_0001,
                                                 exact_mass,
                                                 C,
                                                 H,
                                                 N,
                                                 O,
                                                 P,
                                                 S,
                                                 valence,
                                                 valence_atoms,
                                                 atoms_available,
                                                 bond_types,
                                                 dummies
                                          FROM substructures WHERE valence <= 4""")

        for i, row in enumerate(test_db_cursor.fetchall()):
            if i == 0:
                self.assertEqual(row, ('*N', 1, 3, 16, 16.0187, 16.018724, 0, 2, 1, 0, 0, 0, 1,
                                       '{0: 1}', 1, '{0: [1.0]}', '[1]'))

        self.assertEqual(i, 8)

        test_db.close()

    def test_calculate_hydrogen_rearrangements(self):

        fragment_ions = [(('C', False), ('C', False)), (('C', False), ('C', False), ('C', False), ('C', True)),
                         (('P', False), ('C', False), ('P', False), ('P', True)),
                         (('C', True), ('C', False), ('C', False), ('C', False)),
                         (('C', False), ('P', False), ('P', True), ('P', False)),
                         (('C', False), ('N', False), ('C', False), ('N', False)),
                         (('N', False), ('C', False), ('N', False), ('C', False)),
                         (),
                         (('C', False),)]
        positive_results = [{0, -2}, {1, 3, -5, -3, -1}, {1, 3, 5, -5, -3, -1}, {1, 3, -5, -3, -1},
                            {1, 3, -5, -3, -1}, {0, 2, -4, -2}, {0, 2, 4, -2}, {0}, {-1}]
        negative_results = [{0, 2, -2}, {1, 3, 5, -5, -3, -1}, {1, 3, -5, -3, -1}, {1, 3, 5, -5, -3, -1},
                            {1, 3, 5, -5, -3, -1}, {0, 2, 4, -4, -2}, {0, 2, 4, -2}, {0}, {1, -1}]

        for fragment_ion, positive_result, negative_result in zip(fragment_ions, positive_results, negative_results):
            self.assertEqual(calculate_hydrogen_rearrangements(fragment_ion, "+"), positive_result)
            self.assertEqual(calculate_hydrogen_rearrangements(fragment_ion, "-"), negative_result)


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

        shutil.copytree(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../test_data"),
                        cls.to_test_results("test_data"))

        cls.maxDiff = None

    def test_init(self):

        db = SubstructureDb(self.to_test_data("substructures.sqlite"),
                            self.to_test_data("connectivity.sqlite"))

        db.cursor.execute("SELECT * FROM substructures")
        first_row = db.cursor.fetchone()[0:18]
        self.assertEqual(first_row, (1, '*:c(:*)CCN', 4, 10, 56, 56.05, 56.05002399999998, 3, 6, 1, 0, 0, 0, 2,
                                     '{3: 2}', 1, '{3: [1.5, 1.5]}', '[4, 5]'))

        self.assertTrue(Chem.MolFromSmiles(first_row[1], False))
        self.assertEqual(len(db.cursor.fetchall()), 140)

        db.cursor.execute("SELECT * FROM hmdbid_substructures")
        first_row = db.cursor.fetchone()
        self.assertEqual(first_row, ('HMDB0000073', 1))
        self.assertEqual(len(db.cursor.fetchall()), 150)

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

        self.assertEqual(db.get_single_edge([3, 4, 2]), {2: {2: None, 3: 1, 4: 1}, 3: {3: None, 4: 1}, 4: {4: None}})

        std = db.generate_substructure_network(min_node_weight=2, return_networkx=True)

        db.cursor.execute("SELECT * FROM filtered_hmdbid_substructures")
        for hmdb in db.cursor.fetchall():

            self.assertTrue(hmdb[1] in std.nodes)

        db.cursor.execute("SELECT DISTINCT substructure_id FROM filtered_hmdbid_substructures")
        self.assertEqual(len(db.cursor.fetchall()), 10)
        self.assertEqual(std.number_of_nodes(), 10)

        self.assertEqual(std.number_of_edges(), 24)

        edge_count = []
        db.cursor.execute("SELECT * FROM substructure_graph")
        for edge in db.cursor.fetchall():
            edge_count.append(std.get_edge_data(edge[0], edge[1])["weight"])

        self.assertEqual(sum(edge_count), 48)

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

        self.assertEqual(len(ests), 49)
        self.assertEqual(len(exacts), 65)

        for exact in exacts:
            self.assertTrue(round(exact) in ests)

        self.assertEqual(db.select_mass_values("0_0001", [120, 87, 87], None),
                         [[120.0423], [87.0446], [87.0446]])
        self.assertEqual(db.select_mass_values("0_0001", [55, 80, 107], None),
                         [[55.0184, 55.0422], [80.0262], [107.0497, 107.0735]])

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
        self.assertEqual(len(db.select_substructures([[3, 6, 1, 0, 0, 0]], None)[0]), 1)
        self.assertEqual(list(db.select_substructures([[3, 6, 1, 0, 0, 0]], None)[0][0].keys()),
                         ['smiles', 'mol', 'bond_types', 'degree_atoms', 'valence', 'atoms_available', 'dummies'])

        substructures = list(db.select_substructures([[3, 6, 1, 0, 0, 0]], None)[0][0].values())
        self.assertEqual([item for i, item in enumerate(substructures) if i != 1],
                         ['*:c(:*)CCN', {3: [1.5, 1.5]}, {3: 2}, 2, 1, [4, 5]])

        self.assertEqual(len(db.select_substructures([[2, 2, 0, 2, 0, 0]], None)[0]), 2)
        self.assertEqual(list(db.select_substructures([[2, 2, 0, 2, 0, 0]], None)[0][0].keys()),
                         ['smiles', 'mol', 'bond_types', 'degree_atoms', 'valence', 'atoms_available', 'dummies'])
        substructures = list(db.select_substructures([[2, 2, 0, 2, 0, 0]], None)[0][0].values())
        self.assertEqual([item for i, item in enumerate(substructures) if i != 1],
                         ['*:c(O)c(:*)O', {1: [1.5], 3: [1.5]}, {1: 1, 3: 1}, 2, 2, [0, 5]])

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

    def test_calculate_possible_hydrogenations(self):  # also tests insert_substructure_ion

        records = [self.to_test_data(r + ".xml") for r in ["HMDB0001245", "HMDB0000263"]]

        create_substructure_database(records, self.to_test_results("substructures.sqlite"), 3, 20, method="exhaustive",
                                     isomeric_smiles=True, max_degree=6, max_atoms_available=2)

        db = SubstructureDb(self.to_test_results("substructures.sqlite"))

        search_statement = """SELECT smiles, hydrogen_modification, valence, substructure_ions.substructure_id
                                  FROM substructure_ions 
                                  LEFT JOIN substructures ON substructure_ions.substructure_id = substructures.substructure_id 
                                  WHERE modified_exact_mass__0_0001 > ({} {} 1.007276) - 0.01 
                                  AND modified_exact_mass__0_0001 < ({} {} 1.007276) + 0.01 
                                  AND ion_mode_positive = {}
                           """

        # HMDB0001245 - 2'-deoxycytidine 5'-diphosphate
        mzs = [256.9616, 158.9248, 96.9691, 78.9585]
        h_mods = [-2, -1, 1, -1]
        valences = [2, 1, 1, 1]
        smiles = ["*[C@H]1C[C@H](O)[C@@H](COP(=O)(O)OP(*)(=O)O)O1", "*P(=O)(O)OP(=O)(O)O", "*OP(=O)(O)O", "*P(=O)(O)O"]

        for mz, h_mod, valence, smile in zip(mzs, h_mods, valences, smiles):

            db.cursor.execute(search_statement.format(mz, "+", mz, "+", 0))
            substructure_found = False

            for substructure in db.cursor.fetchall():

                if substructure[0] == smile:

                    substructure_found = True

                    self.assertEqual(substructure[1], h_mod)
                    self.assertEqual(substructure[2], valence)

            self.assertTrue(substructure_found)

        # HMDB0000263 - Phospho(enol)pyruvic acid
        mzs = [62.9628, 64.9785, 80.9734, 89.0233, 94.9892, 98.9841, 104.9735, 116.973611, 122.9842, 140.9947, 150.9791]
        h_mods = [-1, 1, -1, 1, 2, 1, -2, -1, -1, 3, -1]
        valences = [3, 3, 1, 1, 4, 1, 2, 3, 1, 3, 1]
        smiles = ['*OP(*)(*)=O', '*OP(*)(*)=O', "*P(=O)(O)O", "*OC(=C)C(=O)O", "*C(=*)OP(*)(=O)O", "*OP(=O)(O)O",
                  "*C(=C)OP(*)(=O)O", "*C(=O)C(=C)OP(*)(*)=O", "*C(=C)OP(=O)(O)O",
                  "*P(=O)(O)OC(=*)C(=O)O", "*P(=O)(O)OC(=C)C(=O)O"]

        for mz, h_mod, valence, smile in zip(mzs, h_mods, valences, smiles):

            db.cursor.execute(search_statement.format(mz, "-", mz, "-", 1))
            substructure_found = False

            for substructure in db.cursor.fetchall():

                if substructure[0] == smile:

                    substructure_found = True

                    self.assertEqual(substructure[1], h_mod)
                    self.assertEqual(substructure[2], valence)

                    break

            self.assertTrue(substructure_found)


if __name__ == '__main__':
    unittest.main()
