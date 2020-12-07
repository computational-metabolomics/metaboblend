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
import sys
import unittest
import shutil
import tempfile
from metaboblend.databases import *


class IsomorphDbTestCase(unittest.TestCase):
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

    def test_create_connectivity_database(self):

        pkg_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

        if sys.platform == "win32" or sys.platform == "win64":  # TODO: add RI as dependency
            self.path_ri = os.path.join(pkg_path, "tools", "RI_win", "RI3.6-release", "ri36")

        else:

            if sys.platform == "darwin":
                self.path_ri = os.path.join(pkg_path, "tools", "RI_mac", "RI3.6-release", "ri36")

            elif sys.platform == "linux2":
                if "bb" in "socket.gethostname":
                    self.path_ri = os.path.join(pkg_path, "tools", "RI_unix", "RI3.6-release", "ri36")
                else:
                    self.path_ri = os.path.join(pkg_path, "tools", "RI_bb", "RI3.6-release", "ri36")

            elif sys.platform == "linux":
                self.path_ri = os.path.join(pkg_path, "tools", "RI_unix", "RI3.6-release", "ri36")

            create_connectivity_database(self.to_test_results("connectivity.sqlite"),
                                         3,  # sizes
                                         [1, 2],  # boxes
                                         self.path_ri
                                         )

            ref_db = sqlite3.connect(self.to_test_data("connectivity.sqlite"))
            ref_db_cursor = ref_db.cursor()
            ref_db_cursor.execute("SELECT * FROM subgraphs")

            test_db = sqlite3.connect(self.to_test_results("connectivity.sqlite"))
            test_db_cursor = test_db.cursor()
            test_db_cursor.execute("SELECT * FROM subgraphs")

            test_rows = {}
            for row in test_db_cursor.fetchall():
                test_rows[row[0]] = row

            for row in ref_db_cursor.fetchall():
                self.assertEqual(row, test_rows[row[0]])

            ref_db.close()
            test_db.close()


if __name__ == '__main__':
    unittest.main()
