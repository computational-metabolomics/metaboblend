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
import zipfile
import pickle
from metaboblend.databases import *


def to_test_result(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_results", *args)


class DatabasesTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        os.mkdir(to_test_result())

        zip_ref = zipfile.ZipFile(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                               "data",
                                               "connectivity.zip"
                                               ), 'r')
        zip_ref.extractall(to_test_result())
        zip_ref.close()

        zip_ref = zipfile.ZipFile(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                               "data",
                                               "test_mols.zip"
                                               ), 'r')
        zip_ref.extractall(to_test_result())
        zip_ref.close()

        zip_ref = zipfile.ZipFile(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                               "data",
                                               "substructures.zip"
                                               ), 'r')
        zip_ref.extractall(to_test_result())
        zip_ref.close()

    def test_reformat_xml(self):
        reformat_xml(to_test_result("test_mols", "HMDB0000073.xml"))

        with open(to_test_result("test_mols", "HMDB0000073.xml"), "r", encoding="utf-8") as fn_hmdb:
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

        with open(to_test_result("test_mols", "parsed_records.dictionary"), "rb") as parsed:
            parsed_records = pickle.load(parsed)

        for hmdb_xml in os.listdir(to_test_result("test_mols", "hmdb")):
            for record_out in parse_xml(to_test_result("test_mols", "hmdb", hmdb_xml)):
                ref = [i for i, hmdb in enumerate(hmdbs) if hmdb == hmdb_xml.replace(".xml", "")][0]

                self.assertEqual(len(record_out), lengths[ref])
                self.assertEqual(record_out["accession"], hmdbs[ref])
                self.assertEqual(record_out["smiles"], smis[ref])
                self.assertEqual(record_out["chemical_formula"], elements[ref])
                self.assertEqual(record_out, parsed_records[hmdb_xml])

    def test_filter_records(self):
        with open(to_test_result("test_mols", "parsed_records.dictionary"), "rb") as p:
            parsed_records = pickle.load(p)

        with open(to_test_result("test_mols", "test_hmdbs.dictionary"), "rb") as test_hmdbs:
            filtered_records = pickle.load(test_hmdbs)

        record_gen = filter_records(parsed_records.values())
        test_filtered_records = {}
        for record in record_gen:
            del record["mol"]
            test_filtered_records[record["HMDB_ID"]] = record

        self.assertEqual(test_filtered_records, filtered_records)

    def test_get_substructure_bond_idx(self):
        with open(to_test_result("test_mols", "test_hmdbs.dictionary"), "rb") as test_hmdbs:
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
        with open(to_test_result("test_mols", "test_hmdbs.dictionary"), "rb") as test_hmdbs:
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
        with open(to_test_result("test_mols", "test_hmdbs.dictionary"), "rb") as test_hmdbs:
            record_dict = pickle.load(test_hmdbs)["HMDB0000186"]
            mol = Chem.MolFromSmiles(record_dict["smiles"])

        libs = [{'smiles': '*[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O', 'bond_types': {4: [1.0]},
                 'degree_atoms': {4: 1}, 'valence': 1, 'atoms_available': 1, 'dummies': [5]},
                {'smiles': '*O[C@@H]1O[C@H](CO)[C@H](*)[C@H](O)[C@H]1O', 'bond_types': {5: [1.0], 11: [1.0]},
                 'degree_atoms': {5: 1, 11: 1}, 'valence': 2, 'atoms_available': 2, 'dummies': [6, 12]}]

        for edges in [[(0, 1, 2, 22, 3, 16, 17, 18, 19, 20, 21), (0, 1, 2, 22, 3, 4, 16, 17, 18, 19, 20)]]:
            for i, edge_idx in enumerate(edges):
                lib = get_substructure(mol, edge_idx)
                del lib["mol"]
                self.assertEqual(lib, libs[i])

    def test_get_elements(self):
        compositions = [{'C': 8, 'H': 11, 'N': 1, 'O': 2, 'P': 0, 'S': 0, '*': 0},
                        {'C': 6, 'H': 12, 'N': 0, 'O': 6, 'P': 0, 'S': 0, '*': 0},
                        {'C': 9, 'H': 11, 'N': 1, 'O': 3, 'P': 0, 'S': 0, '*': 0},
                        {'C': 12, 'H': 22, 'N': 0, 'O': 11, 'P': 0, 'S': 0, '*': 0}]

        with open(to_test_result("test_mols", "test_hmdbs.dictionary"), "rb") as test_hmdbs:
            record_dicts = pickle.load(test_hmdbs)
            for i, record_dict in enumerate(record_dicts.values()):
                mol = Chem.MolFromSmiles(record_dict["smiles"])

                self.assertEqual(get_elements(mol), compositions[i])

    def test_calculate_exact_mass(self):
        masses = [153.07897899999998, 180.06338999999997, 181.07389399999997, 342.1162150000005]

        with open(to_test_result("test_mols", "test_hmdbs.dictionary"), "rb") as test_hmdbs:
            record_dicts = pickle.load(test_hmdbs)
            for i, record_dict in enumerate(record_dicts.values()):
                mol = Chem.MolFromSmiles(record_dict["smiles"])

                self.assertEqual(calculate_exact_mass(mol), masses[i])

        ref_db = sqlite3.connect(to_test_result("substructures.sqlite"))
        ref_db_cursor = ref_db.cursor()
        ref_db_cursor.execute("SELECT exact_mass__0_0001, mol FROM substructures")
        for row in ref_db_cursor.fetchall():
            self.assertEqual(round(calculate_exact_mass(Chem.Mol(row[1])), 4), row[0])

        ref_db.close()

    def test_update_substructure_database(self):  # requires create_compound_database from SubstructureDb
        db = SubstructureDb(to_test_result("test_db.sqlite"), "")
        db.create_compound_database()
        db.close()

        records = os.listdir(to_test_result("test_mols", "hmdb"))
        for record in records:
            update_substructure_database(to_test_result("test_mols", "hmdb", record),
                                         to_test_result("test_db.sqlite"), 3, 7, method="exhaustive")

        test_db = sqlite3.connect(to_test_result("test_db.sqlite"))
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

        self.assertEqual(total_rows, 585)

        test_db_cursor.execute("SELECT * FROM hmdbid_substructures")
        for i, row in enumerate(test_db_cursor.fetchall()):
            if i == 0:
                self.assertEqual(row, ('HMDB0000073', '*:c(:*)CCN'))
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

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(to_test_result()):
            if os.path.isdir(to_test_result("connectivity")):
                if os.path.isdir(to_test_result("connectivity", "pkls")):
                    for pkl in os.listdir(to_test_result("connectivity", "pkls")):
                        if os.path.isfile(to_test_result("connectivity", "pkls", pkl)):
                            os.remove(to_test_result("connectivity", "pkls", pkl))

                    os.rmdir(to_test_result("connectivity", "pkls"))

                if os.path.isfile(to_test_result("connectivity", "k_graphs.sqlite")):
                    os.remove(to_test_result("connectivity", "k_graphs.sqlite"))

                os.rmdir(to_test_result("connectivity"))

            if os.path.isdir(to_test_result("test_mols")):
                if os.path.isfile(to_test_result("test_mols", "test_hmdbs.dictionary")):
                    os.remove(to_test_result("test_mols", "test_hmdbs.dictionary"))

                if os.path.isfile(to_test_result("test_mols", "parsed_records.dictionary")):
                    os.remove(to_test_result("test_mols", "parsed_records.dictionary"))

                if os.path.isfile(to_test_result("test_mols", "HMDB0000073.xml")):
                    os.remove(to_test_result("test_mols", "HMDB0000073.xml"))

                if os.path.isdir(to_test_result("test_mols", "hmdb")):
                    for hmdb_xml in os.listdir(to_test_result("test_mols", "hmdb")):
                        os.remove(to_test_result("test_mols", "hmdb", hmdb_xml))

                    os.rmdir(to_test_result("test_mols", "hmdb"))

                os.rmdir(to_test_result("test_mols"))

            if os.path.isfile(to_test_result("test_db.sqlite")):
                os.remove(to_test_result("test_db.sqlite"))

            if os.path.isfile(to_test_result("substructures.sqlite")):
                os.remove((to_test_result("substructures.sqlite")))

            os.rmdir(to_test_result())


if __name__ == '__main__':
    unittest.main()
