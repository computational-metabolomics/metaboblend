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
import sys
import unittest
import zipfile
from metaboverse import *

sys.path.append(os.path.join("..", "metaboverse"))
from build_structures import *


def to_test_result(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_results", *args)


class BuildStructuresTestCase(unittest.TestCase):

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

    def test_build(self):
        db = SubstructureDb(to_test_result("substructures_copy.sqlite"), to_test_result("connectivity", "pkls"),
                            to_test_result("connectivity", "k_graphs.sqlite"))

    def test_build_from_subsets(self):
        pass

    def test_gen_subs_table(self):
        pass

    def test_subset_sum(self):
        pass

    def test_combine_ecs(self):
        pass

    def test_reindex_atoms(self):
        pass

    def test_add_bonds(self):
        pass

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
