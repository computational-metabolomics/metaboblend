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
from io import BytesIO
import zipfile
import tempfile
import pickle
from metaboblend.auxiliary import *


class AuxiliaryTestCase(unittest.TestCase):
    temp_results_dir = None
    temp_results_name = None

    @classmethod
    def to_test_result(cls, *args):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), cls.temp_results_name, *args)

    @classmethod
    def setUpClass(cls):
        cls.temp_results_dir = tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.realpath(__file__)))
        cls.temp_results_name = cls.temp_results_dir.name

        cls.lines_geng = [b'E?oo', b'ECO_', b'ECQ_', b'ECZ?', b'ECX_', b'ECYO', b'EEh_', b'EQhO']

        for compr_data in ["connectivity.zip", "test_mols.zip", "substructures.zip"]:
            zip_ref = zipfile.ZipFile(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                                   "data",
                                                   compr_data
                                                   ), 'r')
            zip_ref.extractall(cls.to_test_result())
            zip_ref.close()

        with open(cls.to_test_result("test_aux", "mappings.pkl"), "rb") as mappings_pkl:
            cls.mappings = pickle.load(mappings_pkl)

        with open(cls.to_test_result("test_aux", "gi_out.pkl"), "rb") as gi_out_pkl:
            cls.gi_out = pickle.load(gi_out_pkl)

        cls.p_list = []
        for G, p in calculate_complete_multipartite_graphs([1, 2], 3):
            cls.p_list.append(p)
            cls.final_graph = G

    def test_calculate_complete_multipartite_graphs(self):
        self.assertEqual(self.p_list, [(1, 1), (1, 2), (2, 2), (1, 1, 1), (1, 1, 2), (1, 2, 2), (2, 2, 2)])
        self.assertEqual(nx.number_of_edges(self.final_graph), 12)
        self.assertEqual(nx.number_of_nodes(self.final_graph), 6)
        self.assertEqual(list(self.final_graph.edges()), [(0, 2), (0, 3), (0, 4), (0, 5), (1, 2), (1, 3), (1, 4),
                                                          (1, 5), (2, 4), (2, 5), (3, 4), (3, 5)])
        self.assertEqual(list(self.final_graph.nodes()), [0, 1, 2, 3, 4, 5])

    def test_draw_subgraph(self):
        # INSERT: 1 CU 1 2 (2, 2) (3, 3) ((1, 2), (2, 1)) 4 3
        small_plt, small_sG = draw_subgraph([(0, 2), (1, 2), (1, 3)], ((1, 2), (2, 1)))
        # INSERT: 1 DQo 2 3 (1, 2, 2) (2, 4, 2) ((2,), (2, 2), (1, 1)) 5 4
        large_plt, large_sG = draw_subgraph([(0, 1), (0, 2), (1, 3), (2, 4)], ((2,), (2, 2), (1, 1)))

        self.assertEqual(list(small_sG.nodes()), [0, 2, 1, 3])
        self.assertEqual(small_sG.number_of_nodes(), 4)
        self.assertEqual(list(small_sG.edges()), [(0, 2), (2, 1), (1, 3)])
        self.assertEqual(small_sG.number_of_edges(), 3)
        self.assertEqual(list(large_sG.nodes()), [0, 1, 2, 3, 4])
        self.assertEqual(large_sG.number_of_nodes(), 5)
        self.assertEqual(list(large_sG.edges()), [(0, 1), (0, 2), (1, 3), (2, 4)])
        self.assertEqual(large_sG.number_of_edges(), 4)

    def test_graph_to_ri(self):
        k_graph = graph_to_ri(self.final_graph, "k_graph")
        self.assertEqual(nx.number_of_nodes(self.final_graph) + nx.number_of_edges(self.final_graph) + 3,
                         k_graph.count("\n"))

        for i, line_geng in enumerate(self.lines_geng):
            sG = nx.read_graph6(BytesIO(line_geng))
            subgraph = graph_to_ri(sG, "subgraph")

            self.assertEqual(nx.number_of_nodes(sG) + nx.number_of_edges(sG) + 3, subgraph.count("\n"))

    def test_iso_aux(self):  # tests graph_info, valences, sort_subgraphs
        for line_geng, mappings, gi_val in zip(self.lines_geng, self.mappings, self.gi_out):
            sG = nx.read_graph6(BytesIO(line_geng))

            if len(mappings) > 0:
                gi = graph_info(self.p_list[-1], sG, mappings, )
                self.assertEqual(gi, gi_val)  # test gi vs reference values

                for m in mappings:
                    ug = nx.relabel_nodes(sG, m, copy=True)
                    vn = get_degrees(self.p_list[-1], ug)
                    self.assertEqual(len(vn), len(self.p_list[-1]))
                    [self.assertEqual(len(t), 2) for t in vn]

                for vn in gi:
                    sorted_subgraphs = sort_subgraphs(gi[vn])
                    self.assertLessEqual(len(sorted_subgraphs), len(gi[vn]))

                    for subgraph in gi[vn]:
                        self.assertTrue(sorted([tuple(sorted(e)) for e in subgraph]) in sorted_subgraphs)


if __name__ == '__main__':
    unittest.main()
