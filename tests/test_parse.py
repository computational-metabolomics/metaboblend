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
import shutil
import tempfile
from metaboblend.parse import *


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

    def test_parse_msp(self):
        for i, ms in enumerate(parse_msp(self.to_test_data("mona_msp.msp"))):

            if i < 2:
                self.assertEqual(ms, None)
            else:
                print(ms)

        self.assertEqual(list(parse_msp(self.to_test_data("massbank_msp.txt")))[0], None)


if __name__ == '__main__':
    unittest.main()
