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


import unittest
import zipfile
import pickle
from metaboblend.build_structures import *
from metaboblend.databases import *


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
        db = SubstructureDb(to_test_result("substructures.sqlite"), to_test_result("connectivity", "pkls"),
                            to_test_result("connectivity", "k_graphs.sqlite"))
        smis = [{'NCCc1cc(O)ccc1O', 'NCCc1cccc(O)c1O', 'NCCc1cc(O)cc(O)c1', 'NCCc1ccc(O)c(O)c1'},
                None,
                {'N[C@@H](Cc1ccc(O)cc1)C(=O)O', 'N[C@@H](Cc1cccc(O)c1)C(=O)O', 'N[C@H](Cc1ccc(O)cc1)C(=O)O'},
                None]

        std_lens = [15, 76, 5, 1920]

        fragments = [58.005, 60.021, 58.005, 58.005]

        exp_lens = [1, 41, 0, 0]

        with open(to_test_result("test_mols", "test_hmdbs.dictionary"), "rb") as test_hmdbs:
            record_dicts = pickle.load(test_hmdbs)
            for i, record_dict in enumerate(record_dicts.values()):
                record_dict["mol"] = Chem.MolFromSmiles(record_dict["smiles"])

                build(mc=[record_dict["C"], record_dict["H"], record_dict["N"], record_dict["O"], record_dict["P"],
                          record_dict["S"]], exact_mass=record_dict["exact_mass"],
                      fn_out=to_test_result(record_dict["HMDB_ID"] + ".smi"), heavy_atoms=range(4, 9),
                      max_valence=4, accuracy="1", max_atoms_available=2, max_n_substructures=3,
                      path_connectivity_db=to_test_result("connectivity", "k_graphs.sqlite"),
                      path_pkls=to_test_result("connectivity", "pkls"), path_substructure_db=to_test_result("substructures.sqlite"))

                j = 0
                unique_smis = set()
                with open(to_test_result(record_dict["HMDB_ID"] + ".smi"), "r") as smi_out:
                    for line in smi_out:
                        j += 1
                        unique_smis.add(line.split()[0])

                self.assertEqual(j, std_lens[i])

                if smis[i] is not None:
                    self.assertEqual(unique_smis, smis[i])
                else:
                    self.assertTrue(len(unique_smis) == 51 or len(unique_smis) == 1892)

                if os.path.isfile(to_test_result(record_dict["HMDB_ID"] + ".smi")):
                    os.remove(to_test_result(record_dict["HMDB_ID"] + ".smi"))

                build(mc=[record_dict["C"], record_dict["H"], record_dict["N"], record_dict["O"], record_dict["P"],
                          record_dict["S"]], exact_mass=record_dict["exact_mass"],
                      fn_out=to_test_result(record_dict["HMDB_ID"] + ".smi"), heavy_atoms=range(4, 9), max_valence=4,
                      accuracy="1", max_atoms_available=2, max_n_substructures=3, fragment_mass=fragments[i], ppm=15,
                      path_connectivity_db=to_test_result("connectivity", "k_graphs.sqlite"),
                      path_pkls=to_test_result("connectivity", "pkls"), path_substructure_db=to_test_result("substructures.sqlite"))

                j = 0
                unique_smis = set()
                with open(to_test_result(record_dict["HMDB_ID"] + ".smi"), "r") as smi_out:
                    for line in smi_out:
                        j += 1
                        unique_smis.add(line.split()[0])

                if i == 0:
                    self.assertEqual(unique_smis, {'NCCc1ccc(O)c(O)c1'})

                self.assertEqual(len(unique_smis), exp_lens[i])

                if os.path.isfile(to_test_result(record_dict["HMDB_ID"] + ".smi")):
                    os.remove(to_test_result(record_dict["HMDB_ID"] + ".smi"))

        db.close()

    def test_gen_subs_table(self):
        db = SubstructureDb(to_test_result("substructures.sqlite"), "")
        table_name = gen_subs_table(db, range(5, 7), 4, 2, 500)

        i = 0
        db.cursor.execute("SELECT heavy_atoms, valence, atoms_available FROM %s" % table_name)
        for row in db.cursor.fetchall():
            i += 1
            self.assertTrue(row[0] in range(5, 7))
            self.assertTrue(row[1] <= 4)
            self.assertTrue(row[2] <= 2)

        self.assertEqual(i, 58)

        db.close()

    def test_subset_sum(self):
        self.assertEqual([s_sum for s_sum in subset_sum([1, 2, 3, 4], 5)], [[2, 3], [1, 4]])

        self.assertEqual(len(list(subset_sum(list(range(60)), 70, 3))), 378)
        self.assertEqual(len(list(subset_sum(list(range(60)), 70, 1000))), 29884)

    def test_combine_ecs(self):
        db = SubstructureDb(to_test_result("substructures.sqlite"), "")
        self.assertEqual(combine_ecs([54.0106, 69.0578], db, "substructures", "0_0001"),
                         [[(3, 2, 0, 1, 0, 0)], [(4, 7, 1, 0, 0, 0)]])
        self.assertEqual(combine_ecs([54, 69], db, "substructures", "1"),
                         [[(3, 2, 0, 1, 0, 0)], [(4, 7, 1, 0, 0, 0), (4, 5, 0, 1, 0, 0)]])
        self.assertEqual(combine_ecs([54.0101, 69.0580], db, "substructures", "0_0001", ppm=15),
                         [[(3, 2, 0, 1, 0, 0)], [(4, 7, 1, 0, 0, 0)]])
        self.assertEqual(combine_ecs([54.0101, 69.0580], db, "substructures", "0_0001", ppm=1), [])
        self.assertEqual(combine_ecs([54.0101, 69.0580], db, "substructures", "0_0001"), [])

        db.close()

    def test_reindex_atoms(self):
        llls = [[{'smiles': '*C(*)C(=O)O', 'mol': None, 'bond_types': {1: [1.0, 1.0]}, 'degree_atoms': {1: 2},
                  'valence': 2, 'atoms_available': 1, 'dummies': [0, 2]},
                 {'smiles': 'NCCc1c:*:*:cc1', 'mol': None, 'bond_types': {4: [1.5], 6: [1.5], 7: [1.5]},
                  'degree_atoms': {4: 1, 7: 1}, 'valence': 2, 'atoms_available': 2, 'dummies': [5, 6]}],
                [{'smiles': '*[C@@H](O)[C@@H](*)O', 'mol': None, 'bond_types': {1: [1.0], 3: [1.0]},
                  'degree_atoms': {1: 1, 3: 1}, 'valence': 2, 'atoms_available': 2, 'dummies': [0, 5]},
                 {'smiles': 'OC1**[C@@H](O)[C@H](O)[C@H]1O', 'mol': None, 'bond_types': {0: [1.0], 3: [1.0], 4: [1.0]},
                  'degree_atoms': {0: 1, 4: 1}, 'valence': 2, 'atoms_available': 2, 'dummies': [2, 3]}],
                [{'smiles': '*C[C@H](N)C(=O)O', 'mol': None, 'bond_types': {2: [1.0]}, 'degree_atoms': {2: 1},
                  'valence': 1, 'atoms_available': 1, 'dummies': [3]},
                 {'smiles': '*c1ccc(O)cc1', 'mol': None, 'bond_types': {1: [1.0]}, 'degree_atoms': {1: 1},
                  'valence': 1, 'atoms_available': 1, 'dummies': [0]}]]

        reindexed = [["*C(*)C(=O)O.NCCc1c:*:*:cc1", [1, 10, 13], [0, 2, 11, 12],
                      {1: [1.0, 1.0], 10: [1.5], 12: [1.5], 13: [1.5]}],
                     ["*[C@@H](O)[C@@H](*)O.OC1**[C@@H](O)[C@H](O)[C@H]1O", [1, 3, 6, 10], [0, 5, 8, 9],
                      {1: [1.0], 3: [1.0], 6: [1.0], 9: [1.0], 10: [1.0]}],
                     ["*C[C@H](N)C(=O)O.*c1ccc(O)cc1", [2, 8], [3, 7], {2: [1.0], 8: [1.0]}]]

        for lll, reindex in zip(llls, reindexed):
            lll[0]["mol"] = Chem.MolFromSmiles(lll[0]["smiles"], False)
            lll[1]["mol"] = Chem.MolFromSmiles(lll[1]["smiles"], False)

            mol_comb, atoms_available, atoms_to_remove, bond_types, bond_mismatch = reindex_atoms(lll)
            self.assertEqual([Chem.MolToSmiles(mol_comb), atoms_available, atoms_to_remove, bond_types], reindex)

    def test_add_bonds(self):
        mol_comb = [Chem.MolFromSmiles("*C(*)C(=O)O.NCCc1c:*:*:cc1", False),
                    Chem.MolFromSmiles("*[C@@H](O)[C@@H](*)O.OC1**[C@@H](O)[C@H](O)[C@H]1O", False),
                    Chem.MolFromSmiles("*C[C@H](N)C(=O)O.*c1ccc(O)cc1", False)]

        atoms_available = [[1, 10, 13], [1, 3, 6, 10], [2, 8]]

        bond_types = [{1: [1.0, 1.0], 10: [1.5], 12: [1.5], 13: [1.5]},
                       {1: [1.0], 3: [1.0], 6: [1.0], 9: [1.0], 10: [1.0]},
                       {2: [1.0], 8: [1.0]}]

        edges = [((0, 1), (0, 2)), ((0, 2), (1, 3)), ((0, 1),)]

        mol_out = [None, "*[CH]1(O)OC2**[CH](O)(C(O)C2O)[CH]1(*)O", "*C[CH](N)(C(=O)O)c1(*)ccc(O)cc1"]

        for i in range(len(atoms_available)):
            mol_e = add_bonds(mol_comb[i], edges[i], atoms_available[i], bond_types[i])

            if i == 0:
                self.assertTrue(mol_e is None)
            else:
                self.assertEqual(Chem.MolToSmiles(mol_e.GetMol(), False), mol_out[i])


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

            if os.path.isfile(to_test_result("HMDB0000073.smi")):
                os.remove(to_test_result("HMDB0000073.smi"))

            if os.path.isfile(to_test_result("HMDB0000122.smi")):
                os.remove(to_test_result("HMDB0000122.smi"))

            if os.path.isfile(to_test_result("HMDB0000158.smi")):
                os.remove(to_test_result("HMDB0000158.smi"))

            if os.path.isfile(to_test_result("HMDB0000186.smi")):
                os.remove(to_test_result("HMDB0000186.smi"))

            os.rmdir(to_test_result())


if __name__ == '__main__':
    unittest.main()
