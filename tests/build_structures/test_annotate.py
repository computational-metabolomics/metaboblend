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
import tempfile
import unittest
from rdkit import Chem

from metaboblend.build_structures.annotate import *


class AnnotateTestCase(unittest.TestCase):
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

    def test_generate_structures(self):  # tests vs build

        db = SubstructureDb(self.to_test_data("substructures.sqlite"), self.to_test_data("connectivity.sqlite"))

        fragments = [56.05, 60.0211, 68.0262, 56.0262]

        with open(self.to_test_data("test_hmdbs.dictionary"), "rb") as test_hmdbs:
            record_dicts = pickle.load(test_hmdbs)
            for i, record_dict in enumerate(record_dicts.values()):
                ms_data = {record_dict["HMDB_ID"]: {"mf": [record_dict["C"], record_dict["H"], record_dict["N"],
                                                           record_dict["O"], record_dict["P"], record_dict["S"]],
                                                    "exact_mass": record_dict["exact_mass"]}}

                # test standard building
                returned_smis = list(
                    generate_structures(ms_data, path_substructure_db=self.to_test_data("substructures.sqlite"),
                                        write_csv_output=True, path_out=self.to_test_results(),
                                        max_degree=6, max_atoms_available=2, max_n_substructures=3,
                                        path_connectivity_db=self.to_test_data("connectivity.sqlite"),
                                        minimum_frequency=None, yield_smis=True, isomeric_smiles=True,
                                        retain_substructures=True))

                returned_smis = returned_smis[0][record_dict["HMDB_ID"]]

                build_smis = build(
                    mf=[record_dict["C"], record_dict["H"], record_dict["N"],
                        record_dict["O"], record_dict["P"], record_dict["S"]],
                    exact_mass=record_dict["exact_mass"],
                    max_n_substructures=3, ppm=None, ncpus=None, table_name=None, isomeric_smiles=True,
                    db=db, tolerance=0.0001, prescribed_substructures=None, max_bde=None
                )

                self.assertEqual(set(build_smis.keys()), set(returned_smis))

                ms_data = {record_dict["HMDB_ID"]: {"mf": [record_dict["C"], record_dict["H"], record_dict["N"],
                                                           record_dict["O"], record_dict["P"], record_dict["S"]],
                                                    "exact_mass": record_dict["exact_mass"],
                                                    "prescribed_mass": fragments[i]}}

                # test prescribed building
                returned_smis = list(
                    generate_structures(ms_data, path_substructure_db=self.to_test_data("substructures.sqlite"),
                                        write_csv_output=True, path_out=self.to_test_results(),
                                        max_degree=6, max_atoms_available=2, max_n_substructures=3,
                                        path_connectivity_db=self.to_test_data("connectivity.sqlite"),
                                        minimum_frequency=None, yield_smis=True, isomeric_smiles=True,
                                        retain_substructures=False))

                returned_smis = returned_smis[0][record_dict["HMDB_ID"]]

                prescribed_substructures = get_possible_fragment_ions(fragments[i], db)

                build_smis = build(
                    mf=[record_dict["C"], record_dict["H"], record_dict["N"],
                        record_dict["O"], record_dict["P"], record_dict["S"]],
                    exact_mass=record_dict["exact_mass"], max_n_substructures=3, ppm=0,
                    ncpus=None, table_name=None, isomeric_smiles=True, db=db, tolerance=0.0001,
                    prescribed_substructures=prescribed_substructures, max_bde=None
                )

                self.assertEqual(set(build_smis.keys()), set(returned_smis))

            ms_data = {}
            for i, record_dict in enumerate(record_dicts.values()):
                record_dict["mol"] = Chem.MolFromSmiles(record_dict["smiles"])
                ms_data[record_dict["HMDB_ID"]] = {"mf": [record_dict["C"], record_dict["H"], record_dict["N"],
                                                          record_dict["O"], record_dict["P"], record_dict["S"]],
                                                   "exact_mass": record_dict["exact_mass"],
                                                   "prescribed_masses": None}

            # test building with multiple inputs
            returned_smi_list = list(
                generate_structures(ms_data, path_substructure_db=self.to_test_data("substructures.sqlite"),
                                    write_csv_output=True, path_out=self.to_test_results(),
                                    max_degree=6, max_atoms_available=2, max_n_substructures=3,
                                    path_connectivity_db=self.to_test_data("connectivity.sqlite"),
                                    minimum_frequency=None, yield_smis=True, isomeric_smiles=True,
                                    retain_substructures=False))

            for i, record_dict in enumerate(record_dicts.values()):
                build_smis = build(
                    mf=[record_dict["C"], record_dict["H"], record_dict["N"],
                        record_dict["O"], record_dict["P"], record_dict["S"]],
                    exact_mass=record_dict["exact_mass"],
                    max_n_substructures=3, ppm=None, ncpus=None, table_name=None, isomeric_smiles=True,
                    prescribed_substructures=None, db=db, tolerance=0.0001, max_bde=None
                )

                self.assertEqual(set(build_smis.keys()), set(returned_smi_list[i][record_dict["HMDB_ID"]]))

        db.close()

    def test_annotate_msn(self):  # tests vs build

        db = SubstructureDb(self.to_test_data("substructures.sqlite"))

        overall_lens = [3, 41, 2, 0]
        smis = [{'NCCc1ccc(O)c(O)c1', 'NCCc1cc(O)ccc1O', 'NCCc1cc(O)cc(O)c1'},
                None,
                {'N[C@@H](Cc1cccc(O)c1)C(=O)O', 'N[C@@H](Cc1ccc(O)cc1)C(=O)O'},
                set()]
        freqs = [1, 0, 0, 0]

        fragments = [56.05, 60.0211, 68.0262, 56.0262]

        with open(self.to_test_data("test_hmdbs.dictionary"), "rb") as test_hmdbs:
            record_dicts = pickle.load(test_hmdbs)
            for i, record_dict in enumerate(record_dicts.values()):

                if not os.path.exists(self.to_test_results("annotate")):
                    os.mkdir(self.to_test_results("annotate"))

                ms_data = {record_dict["HMDB_ID"]: {"mf": [record_dict["C"], record_dict["H"], record_dict["N"],
                                                           record_dict["O"], record_dict["P"], record_dict["S"]],
                                                    "exact_mass": record_dict["exact_mass"],
                                                    "neutral_fragment_masses": fragments}}

                # test standard building
                returned_smis = list(annotate_msn(
                    ms_data, max_degree=6, max_atoms_available=2, max_n_substructures=3,
                    write_csv_output=True, retain_substructures=False, path_out=self.to_test_results(),
                    path_connectivity_db=self.to_test_data("connectivity.sqlite"),
                    path_substructure_db=self.to_test_data("substructures.sqlite"),
                    minimum_frequency=None, yield_smis=True, isomeric_smiles=True
                ))

                returned_smis = returned_smis[0][record_dict["HMDB_ID"]]

                self.assertEqual(len([t[1] for t in returned_smis if t[1] > 1]), freqs[i])

                if smis[i] is not None:
                    self.assertEqual(set(t[0] for t in returned_smis), smis[i])

                if i == 0:
                    self.assertEqual(returned_smis[0][1], 3)

            ms_data = {}
            for i, record_dict in enumerate(record_dicts.values()):
                record_dict["mol"] = Chem.MolFromSmiles(record_dict["smiles"])
                ms_data[record_dict["HMDB_ID"]] = {"mf": [record_dict["C"], record_dict["H"], record_dict["N"],
                                                          record_dict["O"], record_dict["P"], record_dict["S"]],
                                                   "exact_mass": record_dict["exact_mass"],
                                                   "neutral_fragment_masses": fragments}

            os.mkdir(self.to_test_results("annotate_multi"))

            # test building with multiple inputs
            returned_smi_list = list(annotate_msn(
                ms_data, max_degree=6, max_atoms_available=2, max_n_substructures=3,
                path_out=self.to_test_results("annotate_multi"), write_csv_output=True,
                path_connectivity_db=self.to_test_data("connectivity.sqlite"),
                path_substructure_db=self.to_test_data("substructures.sqlite"),
                minimum_frequency=None, yield_smis=True,
                isomeric_smiles=True, retain_substructures=False
            ))

            for i, record_dict in enumerate(record_dicts.values()):
                self.assertEqual(len(returned_smi_list[i][record_dict["HMDB_ID"]]), overall_lens[i])

        db.close()


if __name__ == '__main__':
    unittest.main()
