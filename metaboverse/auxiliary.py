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

import itertools
import networkx as nx
import pylab as plt

def calculate_complete_multipartite_graphs(sizes, boxes):
    for nb in range(2, boxes+1):
        for p in itertools.combinations_with_replacement(sizes, nb):
            g = nx.complete_multipartite_graph(*p)
            yield g, p


def cols_dict(sizes):
    num = 0
    temp = [0]

    for s in sizes:
        num += s
        temp.append(num)

    colsDict = {}
    cols = ["b", "r", "y", "g"]

    for i in range(len(temp)-1):
        colsDict[tuple(range(temp[i], temp[i+1]))] = cols[i]

    return colsDict


def valences(sizes, ug):
    num = 0
    temp = [0]
    for s in sizes:
        num += s
        temp.append(num)

    degrees_nodes = ()
    for i in range(len(temp)-1):
        degree_nodes_sub = [t[1] for t in ug.degree(range(temp[i], temp[i + 1]))]
        degrees_nodes += (tuple(degree_nodes_sub),)

    return degrees_nodes


def _hashGraph(g): return tuple(sorted(g.edges()))


def kp_complete_graphs(sizes, boxes):
    for nboxes in range(2, boxes+1):
        for p in itertools.combinations_with_replacement(sizes, nboxes):
            nodes, edges = complete_k_partite(p)
            yield nodes, edges, p


def match_subgraphs(ug_comp, usg, sgi, sizes):
    frags = {}

    for sg in sgi:
        ug_comp.set_vertex_filter(None)
        ug_comp.set_edge_filter(None)

        vmask, emask = mark_subgraph(ug_comp, usg, sg)
        ug_comp.set_filters(emask, vmask)

        vt, vn = valences(sizes, ug_comp)

        edges = [(int(e.target()), int(e.source())) for e in ug_comp.edges()]
        edges.sort()

        if str(vt) not in frags:
            frags[str(vt)] = [edges]

        else:
            if edges not in frags[str(vt)]:
                frags[str(vt)].append(edges)

    return frags


def draw_subgraph(edges, vn):
    plt.title(str(vn))

    sG = nx.Graph()
    sG.add_edges_from([(e[0], e[1]) for e in edges])

    pos = nx.circular_layout(sG)
    nx.draw(sG, pos)

    cols = ["b", "r", "g", "y"]
    cD = {}

    i = 0
    for j, substructure in enumerate(vn):
        if len(substructure) == 1:
            cD[(i,)] = cols[j]
            i += 1

        elif len(substructure) == 2:
            cD[(i, i + 1)] = cols[j]
            i += 2

    for k in cD.keys():
        nx.draw_networkx_nodes(sG, pos=pos, nodelist=k, node_color=cD[k], node_size=800, alpha=1.0)

    nx.draw_networkx_labels(sG, pos=pos)

    plt.show()


def graph_to_ri(graph, name):
    out = "#{}\n".format(name)
    out += "{}\n".format(graph.number_of_nodes())

    for n in graph.nodes():
        out += "n\n"

    out += "{}\n".format(graph.number_of_edges())

    for e in graph.edges():
        out += "{} {} e\n".format(e[0], e[1])

    return out


def graph_info(sizes, sG, mappings):
    frags = {}

    for m in mappings:
        ug = nx.relabel_nodes(sG, m, copy=True)
        vn = valences(sizes, ug)

        e = list(ug.edges())
        e.sort()

        if str(vn) not in frags:
            frags[str(vn)] = [e]
        else:
            if e not in frags[str(vn)]:
                frags[str(vn)].append(e)

    return frags, (sizes, sG, mappings)


def sort_subgraphs(subgraphs):
    sorted_subgraphs = set()

    for fr in subgraphs:
        sorted_subgraphs.add(tuple(sorted([tuple(sorted(e)) for e in fr])))

    return [list(fr) for fr in sorted_subgraphs]
