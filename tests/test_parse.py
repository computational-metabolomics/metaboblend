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


import os
import copy
import shutil
import tempfile
import unittest

from metaboblend.parse import *


class ParseTestCase(unittest.TestCase):
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
        
        cls.neutral_fragment_masses = [155.00332400000002, 173.01262400000002, 175.004724,
                                       250.052324, 251.054324, 252.049224]
        cls.exact_mass = 250.052424
        cls.mf = [10, 10, 4, 2, 0, 1]
        cls.precursor_mz = 251.0597
        cls.fragment_mzs = [156.0106, 174.0199, 176.012, 251.0596, 252.0616, 253.0565]

        cls.maxDiff = None

    def test_parse_msp(self):
        for i, ms in enumerate(parse_msp(self.to_test_data("mona_msp.msp"))):

            if i < 2:
                self.assertEqual(ms, None)
            else:
                self.assertNotEqual(ms, None)

        self.assertEqual(ms, {"ms_id": "AU101101", "mf": self.mf, "precursor_mz": self.precursor_mz,
                              "fragment_mzs": self.fragment_mzs, "precursor_type": "[M+H]+", 'ion_mode': '+',
                              "exact_mass": self.exact_mass, "neutral_fragment_masses": self.neutral_fragment_masses})

        self.assertEqual(list(parse_msp(self.to_test_data("massbank_msp.txt")))[0], None)

        # ensure that parse_msp provides same output as parse_ms_data when providing an msn file
        for parse_msp_dict, parse_ms_dict in zip(parse_msp(self.to_test_data("mona_msp.msp")),
                                                 parse_ms_data(self.to_test_data("mona_msp.msp"))):

            self.assertEqual(parse_msp_dict, parse_ms_dict)

    def test_parse_ms_data(self):

        # exact mass and neutral fragment masses should not be overwritten by parse_ms_data
        full_ms_dict = {"ms_id": "AU101101", "mf": self.mf, "precursor_mz": self.precursor_mz,
                        "fragment_mzs": self.fragment_mzs, "precursor_type": "[M+H]+", "exact_mass": "abcd",
                        "neutral_fragment_masses": ["a", "b", "c", "d"]}

        self.assertEqual(list(parse_ms_data({"AU101101": copy.deepcopy(full_ms_dict)}))[0], full_ms_dict)

        # if exact mass is present should not be overwritten by parse_ms_data
        exact_mass_ms_dict = {"ms_id": "AU101101", "mf": self.mf, "precursor_mz": self.precursor_mz,
                              "fragment_mzs": self.fragment_mzs, "precursor_type": "[M+H]+", "exact_mass": "abc"}

        parsed_exact_mass_ms_dict = list(parse_ms_data({"test": copy.deepcopy(exact_mass_ms_dict)}))[0]
        exact_mass_ms_dict["ms_id"] = "test"
        exact_mass_ms_dict["neutral_fragment_masses"] = self.neutral_fragment_masses
        self.assertEqual(parsed_exact_mass_ms_dict, exact_mass_ms_dict)

        # neutral fragment masses should not be overwritten by parse_ms_data
        neutral_fragment_masses_ms_dict = {"ms_id": "AU101101", "mf": self.mf, "precursor_mz": self.precursor_mz,
                                           "precursor_type": "[M+H]+", "fragment_mzs": self.fragment_mzs,
                                           "neutral_fragment_masses": ["a", "b", "c", "d"]}

        parsed_neutral_fragment_masses_ms_dict = list(parse_ms_data({"AU101101": copy.deepcopy(neutral_fragment_masses_ms_dict)}))[0]
        neutral_fragment_masses_ms_dict["exact_mass"] = self.exact_mass
        neutral_fragment_masses_ms_dict['ion_mode'] = "+"
        self.assertEqual(parsed_neutral_fragment_masses_ms_dict, neutral_fragment_masses_ms_dict)

        uncalculated_ms_dict = {"ms_id": "AU101101", "mf": self.mf, "precursor_mz": self.precursor_mz,
                                "fragment_mzs": self.fragment_mzs, "precursor_type": "[M+H]+"}
        parsed_uncalculated_ms_dict = list(parse_ms_data({"AU101101": copy.deepcopy(uncalculated_ms_dict)}))[0]
        uncalculated_ms_dict["exact_mass"] = self.exact_mass
        uncalculated_ms_dict["neutral_fragment_masses"] = self.neutral_fragment_masses
        uncalculated_ms_dict['ion_mode'] = "+"
        self.assertEqual(parsed_uncalculated_ms_dict, uncalculated_ms_dict)

        # test with msn=False
        generate_structures_dict = {"ms_id": "AU101101", "mf": self.mf, "precursor_mz": self.precursor_mz, 
                                    "prescribed_mass": "m", "precursor_type": "[M+H]+"}
        parsed_generate_structures_dict = list(parse_ms_data({"AU101101": copy.deepcopy(generate_structures_dict)}, False))[0]
        generate_structures_dict["exact_mass"] = self.exact_mass
        generate_structures_dict['ion_mode'] = "+"
        self.assertEqual(parsed_generate_structures_dict, generate_structures_dict)

        # test with exact mass provided
        generate_structures_dict["exact_mass"] = "a"
        parsed_generate_structures_dict = list(parse_ms_data({"AU101101": copy.deepcopy(generate_structures_dict)}, False))[0]
        self.assertEqual(parsed_generate_structures_dict, generate_structures_dict)

    def test_precursor_ions_to_neutral_masses(self):

        ms_dict = {"ms_id": "AU101101", "mf": self.mf, "precursor_mz": self.precursor_mz,
                        "fragment_mzs": self.fragment_mzs, "precursor_type": "[M+H]+"}

        for which in ["both", "fragments", "precursor", "none"]:
            processed_ms_dict = precursor_ions_to_neutral_masses(copy.deepcopy(ms_dict), which)

            if which in ["both", "fragments"]:
                self.assertEqual(processed_ms_dict["neutral_fragment_masses"], self.neutral_fragment_masses)

            if which in ["both", "precursor"]:
                self.assertEqual(processed_ms_dict["exact_mass"], self.exact_mass)

        ms_dict["precursor_type"] = "[M-H]-"

        for which in ["both", "fragments", "precursor", "none"]:
            processed_ms_dict = precursor_ions_to_neutral_masses(copy.deepcopy(ms_dict), which)

            if which in ["both", "fragments"]:
                neutral_fragment_masses = [nfm + 1.007276 for nfm in self.fragment_mzs]
                self.assertEqual(processed_ms_dict["neutral_fragment_masses"], neutral_fragment_masses)

            if which in ["both", "precursor"]:
                self.assertEqual(processed_ms_dict["exact_mass"], self.precursor_mz + 1.007276)

    def test_reformat_msp_input(self):

        unformatted_msp_dict = {'ms_id': 'AU101101', 'mf': 'C10H10N4O2S', 'precursor_mz': '251.0597',
                                'fragment_mzs': self.fragment_mzs,
                                'precursor_type': '[M+H]+'}

        formatted_msp_dict = {'ms_id': 'AU101101', 'mf': self.mf, 'precursor_mz': self.precursor_mz,
                                'fragment_mzs': self.fragment_mzs, 'precursor_type': '[M+H]+',
                                'exact_mass': self.exact_mass, 'ion_mode': '+',
                                'neutral_fragment_masses': self.neutral_fragment_masses}

        self.assertEqual(reformat_msp_input(unformatted_msp_dict), formatted_msp_dict)

        unformatted_msp_dict["precursor_mz"] = None

        with self.assertWarns(UserWarning):
            reformat_msp_input(unformatted_msp_dict)

        unformatted_msp_dict["precursor_mz"] = self.precursor_mz
        unformatted_msp_dict["fragment_mzs"] = []

        with self.assertWarns(UserWarning):
            reformat_msp_input(unformatted_msp_dict)

    def test_mc_to_list(self):

        mc_lists = [[12, 14, 4, 4, 0, 1], [10, 10, 4, 2, 0, 1], [46, 94, 1, 8, 1, 0], [46, 94, 1, 8, 1, 0], None]

        for i, word_formula in enumerate(["C12H14N4O4S", "C10H10N4O2S", "C46H94NO8P", "C46H94NO8P1", "C10H9ClN4O2S"]):
            self.assertEqual(mc_to_list(word_formula), mc_lists[i])


if __name__ == "__main__":
    unittest.main()
