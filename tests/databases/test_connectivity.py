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

import sys
import shutil
import unittest
from io import BytesIO

from metaboblend.databases.connectivity import *


class ConnectivityTestCase(unittest.TestCase):
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

        cls.lines_geng = [b'E?oo', b'ECO_', b'ECQ_', b'ECZ?', b'ECX_', b'ECYO', b'EEh_', b'EQhO']

        with open(cls.to_test_data("mappings.pkl"), "rb") as mappings_pkl:
            cls.mappings = pickle.load(mappings_pkl)

        with open(cls.to_test_data("gi_out.pkl"), "rb") as gi_out_pkl:
            cls.gi_out = pickle.load(gi_out_pkl)

        cls.p_list = []
        for G, p in calculate_complete_multipartite_graphs([1, 2], 3):
            cls.p_list.append(p)
            cls.final_graph = G

        cls.maxDiff = None

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

    def test_create_connectivity_database(self):

        if sys.platform == "linux" or sys.platform == "linux2":

            create_connectivity_database(self.to_test_results("connectivity.sqlite"), 3, [1, 2])

            ref_db = sqlite3.connect(self.to_test_data("connectivity.sqlite"))
            ref_db_cursor = ref_db.cursor()
            ref_db_cursor.execute("SELECT * FROM subgraphs")

            test_db = sqlite3.connect(self.to_test_results("connectivity.sqlite"))
            test_db_cursor = test_db.cursor()
            test_db_cursor.execute("SELECT * FROM subgraphs")

            test_rows = {}
            for row in test_db_cursor.fetchall():
                test_rows[row[0]] = row

            # compare to ref database
            for row in ref_db_cursor.fetchall():

                # check generated graphs are the same
                self.assertEqual(pickle.loads(test_rows[row[0]][9]), pickle.loads(row[9]))

                # check specifications are the same
                self.assertEqual(test_rows[row[0]][0:9], row[0:9])

            ref_db.close()
            test_db.close()


if __name__ == '__main__':
    unittest.main()
