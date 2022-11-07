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
import pickle
import unittest
import tempfile
from rdkit import Chem

from metaboblend.build_structures.annotate import annotate_msn
from metaboblend.databases.results import *


class ResultsDbTestCase(unittest.TestCase):
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

    def test_results_db(self):  # TODO: directly test each unit of ResultsDb

        fragments = [56.05, 60.0211, 68.0262, 56.0262]

        with open(self.to_test_data("test_hmdbs.dictionary"), "rb") as test_hmdbs:
            record_dicts = pickle.load(test_hmdbs)

        ms_data = {}
        for i, record_dict in enumerate(record_dicts.values()):
            record_dict["mol"] = Chem.MolFromSmiles(record_dict["smiles"])
            ms_data[record_dict["HMDB_ID"]] = {"mf": [record_dict["C"], record_dict["H"], record_dict["N"],
                                                      record_dict["O"], record_dict["P"], record_dict["S"]],
                                               "exact_mass": record_dict["exact_mass"],
                                               "neutral_fragment_masses": fragments}

        os.mkdir(self.to_test_results("test_results_db"))
        list(annotate_msn(
            ms_data, max_degree=6, max_atoms_available=2, max_n_substructures=3,
            path_out=self.to_test_results("test_results_db"), write_csv_output=True,
            path_connectivity_db=self.to_test_data("connectivity.sqlite"),
            path_substructure_db=self.to_test_data("substructures.sqlite"),
            minimum_frequency=None, yield_smis=True,
            isomeric_smiles=True, retain_substructures=True
        ))

        # is the sqlite database the size we expect?
        self.assertEqual(os.path.getsize(self.to_test_results("test_results_db", "metaboblend_results.sqlite")), 86016)

        # are the csv files the same as the reference?
        with open(self.to_test_results("test_results_db", "metaboblend_queries.csv"), "r") as results_file, \
                open(self.to_test_data("metaboblend_queries.csv"), "r") as test_file:

            for results_line, test_line in zip(results_file, test_file):
                self.assertEqual(results_line, test_line)

        with open(self.to_test_results("test_results_db", "metaboblend_structures.csv"), "r") as results_file, \
                open(self.to_test_data("metaboblend_structures.csv"), "r") as test_file:

            for results_line, test_line in zip(results_file, test_file):
                self.assertEqual(results_line, test_line)


if __name__ == '__main__':
    unittest.main()
