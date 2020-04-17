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


import unittest
import zipfile
from metaboverse.databases import *


def to_test_result(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_results", *args)


class IsomorphDbTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        os.mkdir(to_test_result())

        if sys.platform == "win32" or sys.platform == "win64":
            cls.path_geng = os.path.join("..", "tools", "nauty25r9_win", "geng")
            cls.path_ri = os.path.join("..", "tools", "RI_win", "RI3.6-release", "ri36")

        elif sys.platform == "darwin":
            cls.path_geng = os.path.join("..", "tools", "nauty25r9_mac", "geng")
            cls.path_ri = os.path.join("..", "tools", "RI_mac", "RI3.6-release", "ri36")

        elif sys.platform == "linux2":
            if "bb" in "socket.gethostname":
                cls.path_geng = os.path.join("..", "tools", "nauty25r9_unix", "geng")
                cls.path_ri = os.path.join("..", "tools", "RI_unix", "RI3.6-release", "ri36")
            else:
                cls.path_geng = os.path.join("..", "tools", "nauty25r9_bb", "geng")
                cls.path_ri = os.path.join("..", "tools", "RI_bb", "RI3.6-release", "ri36")

        elif sys.platform == "linux":
            cls.path_geng = os.path.join("..", "tools", "nauty25r9_unix", "geng")
            cls.path_ri = os.path.join("..", "tools", "RI_unix", "RI3.6-release", "ri36")

        zip_ref = zipfile.ZipFile(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                               "data",
                                               "connectivity.zip"
                                               ), 'r')
        zip_ref.extractall(to_test_result())
        zip_ref.close()

        os.mkdir(to_test_result("pkls"))
        create_isomorphism_database(to_test_result("k_graphs.sqlite"),
                                    to_test_result("pkls"),
                                    3,  # sizes
                                    [1, 2],  # boxes
                                    cls.path_geng,
                                    cls.path_ri
                                    )

    def test_create_isomorphism_database(self):
        ref_db = sqlite3.connect(to_test_result("connectivity", "k_graphs.sqlite"))
        ref_db_cursor = ref_db.cursor()
        ref_db_cursor.execute("SELECT * FROM subgraphs")

        test_db = sqlite3.connect(to_test_result("k_graphs.sqlite"))
        test_db_cursor = test_db.cursor()
        test_db_cursor.execute("SELECT * FROM subgraphs")

        self.assertEqual(ref_db_cursor.fetchall(), test_db_cursor.fetchall())

        ref_db.close()
        test_db.close()

    def test_create_pkls(self):
        for pkl_path in os.listdir(to_test_result("pkls")):
            with open(to_test_result("connectivity", "pkls", pkl_path), "rb") as ref_pkl, \
                 open(to_test_result("pkls", pkl_path), "rb") as test_pkl:

                self.assertEqual(pickle.load(ref_pkl), pickle.load(test_pkl))

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(to_test_result()):
            if os.path.isdir(to_test_result("pkls")):
                for pkl in os.listdir(to_test_result("pkls")):
                    if os.path.isfile(to_test_result("pkls", pkl)):
                        os.remove(to_test_result("pkls", pkl))

                os.rmdir(to_test_result("pkls"))

            if os.path.isfile(to_test_result("k_graphs.sqlite")):
                os.remove(to_test_result("k_graphs.sqlite"))

            if os.path.isdir(to_test_result("connectivity")):
                if os.path.isdir(to_test_result("connectivity", "pkls")):
                    for pkl in os.listdir(to_test_result("connectivity", "pkls")):
                        if os.path.isfile(to_test_result("connectivity", "pkls", pkl)):
                            os.remove(to_test_result("connectivity", "pkls", pkl))

                    os.rmdir(to_test_result("connectivity", "pkls"))

                if os.path.isfile(to_test_result("connectivity", "k_graphs.sqlite")):
                    os.remove(to_test_result("connectivity", "k_graphs.sqlite"))

                os.rmdir(to_test_result("connectivity"))

            os.rmdir(to_test_result())


if __name__ == '__main__':
    unittest.main()
