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
import unittest
import zipfile
import pickle
from metaboverse import *


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
                # print(record_out)
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

        record_gen = filter_records(parsed_records.values(), db_type="hmdb")
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

            os.rmdir(to_test_result())


if __name__ == '__main__':
    unittest.main()
