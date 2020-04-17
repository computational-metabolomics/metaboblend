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
from io import BytesIO
import zipfile
import pickle
from metaboverse.auxiliary import *


def to_test_result(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_results", *args)


class AuxiliaryTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        os.mkdir(to_test_result())

        cls.lines_geng = [b'E?oo', b'ECO_', b'ECQ_', b'ECZ?', b'ECX_', b'ECYO', b'EEh_', b'EQhO']

        zip_ref = zipfile.ZipFile(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                               "data", "test_aux.zip"), 'r')
        zip_ref.extractall(to_test_result())
        zip_ref.close()

        with open(to_test_result("test_aux", "mappings.pkl"), "rb") as mappings_pkl:
            cls.mappings = pickle.load(mappings_pkl)

        with open(to_test_result("test_aux", "gi_out.pkl"), "rb") as gi_out_pkl:
            cls.gi_out = pickle.load(gi_out_pkl)

        cls.p_list = []
        for G, p in calculate_complete_multipartite_graphs([1, 2], 3):
            cls.p_list.append(p)
            cls.final_graph = G

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
                self.assertEqual(gi[0], gi_val)  # test gi vs reference values
                self.assertEqual(gi[1][0], self.p_list[-1])
                self.assertEqual(gi[1][1], sG)
                self.assertEqual(gi[1][2], mappings)

                for m in mappings:
                    ug = nx.relabel_nodes(sG, m, copy=True)
                    vn = valences(self.p_list[-1], ug)
                    self.assertEqual(len(vn), len(self.p_list[-1]))
                    [self.assertEqual(len(t), 2) for t in vn]

                for vn in gi[0]:
                    sorted_subgraphs = sort_subgraphs(gi[0][vn])
                    self.assertLessEqual(len(sorted_subgraphs), len(gi[0][vn]))

                    for subgraph in gi[0][vn]:
                        self.assertTrue(sorted([tuple(sorted(e)) for e in subgraph]) in sorted_subgraphs)

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(to_test_result("test_aux")):

            if os.path.isfile(to_test_result("test_aux", "mappings.pkl")):
                os.remove(to_test_result("test_aux", "mappings.pkl"))

            if os.path.isfile(to_test_result("test_aux", "gi_out.pkl")):
                os.remove(to_test_result("test_aux", "gi_out.pkl"))

            os.rmdir(to_test_result("test_aux"))

        if os.path.isdir(to_test_result()):
            os.rmdir(to_test_result())


if __name__ == '__main__':
    unittest.main()
