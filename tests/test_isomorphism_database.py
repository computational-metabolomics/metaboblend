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
import tempfile
from metaboblend.databases import *


class IsomorphDbTestCase(unittest.TestCase):
    temp_results_dir = None
    temp_results_name = None

    @classmethod
    def to_test_result(cls, *args):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), cls.temp_results_name, *args)

    @classmethod
    def setUpClass(cls):
        cls.temp_results_dir = tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.realpath(__file__)))
        cls.temp_results_name = cls.temp_results_dir.name

        pkg_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        if sys.platform == "win32" or sys.platform == "win64":
            cls.path_ri = os.path.join(pkg_path, "tools", "RI_win", "RI3.6-release", "ri36")

        elif sys.platform == "darwin":
            cls.path_ri = os.path.join(pkg_path, "tools", "RI_mac", "RI3.6-release", "ri36")

        elif sys.platform == "linux2":
            if "bb" in "socket.gethostname":
                cls.path_ri = os.path.join(pkg_path, "tools", "RI_unix", "RI3.6-release", "ri36")
            else:
                cls.path_ri = os.path.join(pkg_path, "tools", "RI_bb", "RI3.6-release", "ri36")

        elif sys.platform == "linux":
            cls.path_ri = os.path.join(pkg_path, "tools", "RI_unix", "RI3.6-release", "ri36")

        zip_ref = zipfile.ZipFile(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                               "data",
                                               "connectivity.zip"
                                               ), 'r')
        zip_ref.extractall(cls.to_test_result())
        zip_ref.close()

        os.mkdir(cls.to_test_result("pkls"))
        create_isomorphism_database(cls.to_test_result("k_graphs.sqlite"),
                                    3,  # sizes
                                    [1, 2],  # boxes
                                    cls.path_ri
                                    )

    def test_create_isomorphism_database(self):
        ref_db = sqlite3.connect(self.to_test_result("connectivity", "k_graphs.sqlite"))
        ref_db_cursor = ref_db.cursor()
        ref_db_cursor.execute("SELECT * FROM subgraphs")

        test_db = sqlite3.connect(self.to_test_result("k_graphs.sqlite"))
        test_db_cursor = test_db.cursor()
        test_db_cursor.execute("SELECT * FROM subgraphs")

        ref_rows = {}
        for row in ref_db_cursor.fetchall():
            ref_rows[row[0]] = row

        for row in test_db_cursor.fetchall():
            self.assertEqual(row, ref_rows[row[0]])

        ref_db.close()
        test_db.close()

    def test_create_pkls(self):
        for pkl_path in os.listdir(self.to_test_result("pkls")):
            with open(self.to_test_result("connectivity", "pkls", pkl_path), "rb") as ref_pkl, \
                 open(self.to_test_result("pkls", pkl_path), "rb") as test_pkl:

                self.assertEqual(pickle.load(ref_pkl), pickle.load(test_pkl))

    @classmethod
    def tearDownClass(cls):
        if cls.temp_results_dir is not None:
            cls.temp_results_dir.cleanup()


if __name__ == '__main__':
    unittest.main()
