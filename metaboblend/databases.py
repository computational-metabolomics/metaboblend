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

import io
import sys
import subprocess
import sqlite3
import tempfile
import pickle
from collections import OrderedDict
import xml.etree.ElementTree as ElementTree
import networkx as nx
from rdkit import Chem
from rdkit.Chem import Recap
from rdkit.Chem import BRICS
from .auxiliary import calculate_complete_multipartite_graphs, graph_to_ri, graph_info, sort_subgraphs


def reformat_xml(source, encoding="utf8"):
    """
    Reformats HMDB xml files to be compatible with :py:meth:`metaboblend.databases.parse_xml`; some such files do not
    contain a `<hmdb xmlns="http://www.hmdb.ca">` header.

    :param source: Path to file to be reformatted.

    :param encoding: Encoding of source file.

    :return: Source file destination.
    """

    with io.open(source, "r", encoding=encoding) as xml:
        xml_contents = xml.readlines()
        if "hmdb" in xml_contents[1]:
            return source

        xml_contents.insert(1, "<hmdb xmlns=\"http://www.hmdb.ca\"> \n")

    with io.open(source, "w", encoding=encoding) as xml:
        xml_contents = "".join(xml_contents)
        xml.write(xml_contents)
        xml.write("</hmdb>")

    return source


def parse_xml(source, encoding="utf8", reformat=False):
    """
    Parses the contents of HMDB xml files to to extract information for the generation of substructures.

    :param source: Source file destination.

    :param encoding: Encoding of source file.

    :param reformat: Whether to apply :py:meth:`metaboblend.databases.reformat_xml` to the XML file. Is required for
        XML files recording single metabolites.

        * **True** Add a `<hmdb xmlns="http://www.hmdb.ca">` header to the XML file before parsing.

        * **False** Parse the XML file as it is (recommended if header is present).

    :return: The XML file converted to a dictionary.
    """

    if reformat:
        reformat_xml(source, encoding)

    with io.open(source, "r", encoding=encoding) as inp:
        record_out = OrderedDict()

        inp.readline()
        inp.readline()

        xml_record = ""
        path = []

        for line in inp:
            xml_record += line
            if line == "</metabolite>\n" or line == "</drug>\n":

                if sys.version_info[0] == 3:
                    inp = io.StringIO(xml_record)
                else:
                    inp = io.BytesIO(xml_record.encode('utf-8').strip())

                for event, elem in ElementTree.iterparse(inp, events=("start", "end")):
                    if event == 'end':
                        path.pop()

                    if event == 'start':
                        path.append(elem.tag)
                        if elem.text is not None:
                            if elem.text.replace(" ", "") != "\n":

                                path_elem = ".".join(map(str, path[1:]))
                                if path_elem in record_out:
                                    if type(record_out[path_elem]) != list:
                                        record_out[path_elem] = [record_out[path_elem]]
                                    record_out[path_elem].append(elem.text)
                                else:
                                    record_out[path_elem] = elem.text

                xml_record = ""
                yield record_out
                record_out = OrderedDict()


class SubstructureDb:
    """
    Methods for interacting with the SQLITE3 substructure and connectivity databases. Provides a connection to the
    substructure database and, if provided, the connectivity database.

    :param path_substructure_db: Path to the substructure database.

    :param path_connectivity_db: Path to the connectivity database.

    :param conn: A :py:meth:`sqlite3.connection` to the substructure database; the connectivity database will be attached
        as 'graphs'.

    :param cursor: A :py:meth:`sqlite3.connection.cursor` for the substructure database; the connectivity database will
        be attached as "graphs".
    """

    def __init__(self, path_substructure_db, path_connectivity_db=None):
        """Constructor method"""

        self.path_substructure_db = path_substructure_db
        self.path_connectivity_db = path_connectivity_db

        self.conn = sqlite3.connect(self.path_substructure_db)
        self.cursor = self.conn.cursor()

        if self.path_connectivity_db is not None:
            self.cursor.execute("ATTACH DATABASE '%s' as 'graphs';" % self.path_connectivity_db)

    def select_compounds(self, cpds=[]):
        """
        Select all keys from the compounds table of the substructure database, filtered by HMDB IDs if provided.

        :param cpds: A list of HMDB IDs by which to filter the compounds; applies no filter if not provided.

        :return: Returns the result of the compounds table query via :py:meth:`sqlite3.connection.cursor.fetchall`;
            essentially provides a list containing a list for the values of each row.
        """

        if len(cpds) > 0:
            sql = " WHERE hmdbid in ('%s')" % ("', '".join(map(str, cpds)))
        else:
            sql = ""

        self.cursor.execute("""SELECT DISTINCT hmdbid, exact_mass, formula, C, H, N, O, P, S, smiles 
                               FROM compounds%s""" % sql)

        return self.cursor.fetchall()

    def filter_hmdbid_substructures(self, min_node_weight=2):
        """
        Filters `hmdbid_substructures` table in the substructure database by the frequency of the substructure; allows
        for the removal of non-unique (min_node_weight = 2) or infrequent (min_node_weight > 2) substructures. Generates
        a new table, *unique_hmdbid*, which contains distinct hmdbids and *filtered_hmdbid_substructures* which contains
        the filtered substructures with their frequency.

        :param min_node_weight: Minimal count of the substructure within 'hmdbid_substructures'.
        """

        self.cursor.execute("DROP TABLE IF EXISTS unique_hmdbid")
        self.cursor.execute("DROP TABLE IF EXISTS filtered_hmdbid_substructures")

        self.cursor.execute("CREATE TABLE unique_hmdbid AS SELECT DISTINCT hmdbid FROM compounds")

        self.cursor.execute("""CREATE TABLE filtered_hmdbid_substructures AS
                                   SELECT smiles, COUNT(*) FROM hmdbid_substructures
                                   GROUP BY smiles HAVING COUNT(*) >=%s""" % min_node_weight)

    def generate_substructure_network(self, method="default", min_node_weight=2, remove_isolated=False):
        """
        Generate networks to explore the co-occurence of substructures in the substructure database.

        :param method: Which method to use for the generation of a substructure co-occurence network.

            * **default** Generate a standard substructure network using the most time-efficient method. Node weights
                represent the frequency of substructures in the substructure database whilst edge weights represent the
                co-occurence of these substructures. See
                :py:meth:`metaboblend.databases.default_substructure_network`.

            * **extended** Alternative method for generating substructure networks that generates an intermediate
                metabolite-substructure linkage network; returns the same network as the **default** method. See
                :py:meth:`metaboblend.databases.extended_substructure_network`.

            * **parent_structure_linkage**  Uses the **extended** method and returns the intermediate
                metabolite-substructure linkage network; instead of weighted edges between substructures being present,
                the substructures are connected to the metabolites by which they were generated.

        :param min_node_weight: Minimum frequency of substructures to be included in the network; a minimum node weight
            of 2 means that all non-unique substructures will be present.

        :param remove_isolated: Isolated nodes are substructures that do not co-occur with any others.

            * **True** Remove isolated nodes from the final network.

            * **False** Allow isolated nodes to remain in the final network.

        :return: A :py:meth:`networkx.Graph` object containing the specified network.
        """

        substructure_graph = nx.Graph()
        self.filter_hmdbid_substructures(min_node_weight)

        self.cursor.execute("SELECT * FROM unique_hmdbid")
        unique_hmdb_ids = self.cursor.fetchall()

        self.cursor.execute("SELECT * FROM filtered_hmdbid_substructures")
        # add node for each unique substructure, weighted by count
        for unique_substructure in self.cursor.fetchall():
            substructure_graph.add_node(unique_substructure[0], weight=unique_substructure[1])

        # generate different flavours of network
        if method == "default":
            substructure_graph = self.default_substructure_network(substructure_graph, unique_hmdb_ids)
        elif method == "extended":
            substructure_graph = self.extended_substructure_network(substructure_graph, unique_hmdb_ids,
                                                                    include_parents=False)
        elif method == "parent_structure_linkage":
            substructure_graph = self.extended_substructure_network(substructure_graph, unique_hmdb_ids,
                                                                    include_parents=True)

        if remove_isolated:
            substructure_graph.remove_nodes_from(list(nx.isolates(substructure_graph)))

        return substructure_graph

    def extended_substructure_network(self, substructure_graph, unique_hmdb_ids, include_parents=False):
        """
        Extended method for substructure network generation - ie, an intermediate substructure network is generated
        which involves substructure nodes, weighted by frequency, and metabolite nodes; unweighted edges are found
        between substructures and metabolites, representing the structures from which substructures were generated.
        This intermediate network can be returned, or it may be transformed into the standard substructure network
        generated by :py:meth:`metaboblend.databases.default_substructure_network`.

        :param substructure_graph: A :py:meth:`networkx.Graph` containing substructure nodes, weighted by their
            frequencies in the substructure database. Can be passed by
            :py:meth:`metaboblend.databases.generate_substructure_network`.

        :param unique_hmdb_ids: A list containing distinct HMDB IDs from the substructure database.

        :param include_parents: Which type of substructure network to return.

            * **True** Return the intermediate network containing metabolite nodes.

            * **False** Return the standard substructure network.

        :return: A :py:meth:`networkx.Graph` object containing the specified network.
        """

        # add node for each parent structure
        for unique_hmdb_id in unique_hmdb_ids:
            substructure_graph.add_node(unique_hmdb_id[0])

        # add edge for each linked parent structure and substructure
        self.cursor.execute("""SELECT * FROM hmdbid_substructures WHERE smiles IN 
                                   (SELECT smiles FROM filtered_hmdbid_substructures)""")
        for hmdbid_substructures in self.cursor.fetchall():
            substructure_graph.add_edge(hmdbid_substructures[0], hmdbid_substructures[1])

        if not include_parents:
            # remove parent structures and replace with linked, weighted substructures
            for unique_hmdb_id in unique_hmdb_ids:
                for adj1 in substructure_graph.adj[unique_hmdb_id[0]]:
                    for adj2 in substructure_graph.adj[unique_hmdb_id[0]]:
                        if substructure_graph.has_edge(adj1, adj2):
                            substructure_graph[adj1][adj2]['weight'] += 1
                        else:
                            substructure_graph.add_edge(adj1, adj2, weight=1)
                substructure_graph.remove_node(unique_hmdb_id[0])

            # remove self-loops
            substructure_graph.remove_edges_from(nx.selfloop_edges(substructure_graph))

        return substructure_graph

    def default_substructure_network(self, substructure_graph, unique_hmdb_ids):
        """
        Standard method for generating substructure co-occurence networks.

        :param substructure_graph: A :py:meth:`networkx.Graph` containing substructure nodes, weighted by their
            frequencies in the substructure database. Can be passed by
            :py:meth:`metaboblend.databases.generate_substructure_network`.

        :param unique_hmdb_ids: A list containing distinct HMDB IDs from the substructure database.

        :return: A :py:meth:`networkx.Graph` object containing the specified network.
        """

        # add edges by walking through hmdbid_substructures
        for unique_hmdb_id in unique_hmdb_ids:
            self.cursor.execute("""SELECT * FROM hmdbid_substructures 
                                       WHERE smiles IN (SELECT smiles FROM filtered_hmdbid_substructures) 
                                       AND hmdbid = '%s'""" % unique_hmdb_id)
            nodes = []
            for substructure in self.cursor.fetchall():
                for node in nodes:
                    if substructure_graph.has_edge(substructure[1], node):
                        substructure_graph[substructure[1]][node]['weight'] += 1
                    else:
                        substructure_graph.add_edge(substructure[1], node, weight=1)

                nodes.append(substructure[1])

        return substructure_graph

    def select_mass_values(self, accuracy, masses, table_name):
        """
        Gets mass values from a table of substructures. Can limit results based on mass at integer level. Used by
        subset_sum to filter substructures prior to building.

        :param accuracy: To which decimal places of accuracy results are to be limited to.

            * **1** Integer level
            * **0_0001** Four decimal places

        :param masses: A list of integers to limit the query by. If a non-empty list of masses is given, the query
            is returned as a list of lists for each mass query, as opposed to a single list of masses.

        :param table_name: Name of the substructure table to be queried.

        :return: Sorted list of mass values from the substructure database, filtered by the supplied parameters.
        """

        if type(masses) == list:
            if len(masses) > 0:
                mass_values = []

                for m in masses:
                    self.cursor.execute("""SELECT DISTINCT exact_mass__{} 
                                           FROM {}
                                           WHERE exact_mass__1 = {}
                                        """.format(accuracy, table_name, m))

                    m_values = [record[0] for record in self.cursor.fetchall()]
                    m_values.sort()

                    mass_values.append(m_values)

                return mass_values

        self.cursor.execute("SELECT DISTINCT exact_mass__{} FROM {}".format(accuracy, table_name))
        mass_values = [record[0] for record in self.cursor.fetchall()]
        mass_values.sort()

        return mass_values

    def select_ecs(self, exact_mass, table_name, accuracy):
        """
        Select elemental compositions based on an exact mass; allows for inexact mass searches based on error (ppm).

        :param exact_mass: The exact mass used to query the substructure database.

        :param table_name: Name of the substructure table to be queried.

        :param accuracy: To which decimal places of accuracy results are to be limited to.

            * **1** Integer level
            * **0_0001** Four decimal places

        :return: Returns the elemental compositions associated with the `exact_mass` via
            :py:meth:`sqlite3.connection.cursor.fetchall`; essentially provides a list containing a list for the values
            of each row (C, H, N, O, P, S).
        """

        self.cursor.execute("""SELECT DISTINCT C, H, N, O, P, S
                                   FROM {} 
                                   WHERE exact_mass__{} = {}
                            """.format(table_name, accuracy, str(exact_mass)))

        return self.cursor.fetchall()

    def k_configs(self, fragment_edges_only=False):
        """
        Obtains strings detailing the valences for each substructure in a connectivity graph and the ID of the related
        PKL file. Used to match a set of substructures to the correct set of non-isomorphic graphs in the connectivity
        database.

        :param fragment_edges_only: If true, only include sets of edges that connect to the fragment ion.

        :return: Dictionary containing the valences as keys and PKL IDs as values.
        """

        self.cursor.execute("""SELECT root, nodes_valences
                                   FROM subgraphs""")

        records = self.cursor.fetchall()
        configs = {}

        for record in records:
            configs[str(record[1])] = []
            fragment_atoms = [i for i in range(len(eval(record[1])[0]))]

            for path in self.paths(pickle.loads(record[0])):
                if fragment_edges_only:
                    all_bonds_connect_to_fragment = True

                    for edge in path:
                        if edge[0] not in fragment_atoms and edge[1] not in fragment_atoms:
                            all_bonds_connect_to_fragment = False

                    # check that all bonds connect to fragment ion
                    if all_bonds_connect_to_fragment:
                        configs[str(record[1])].append(path)

                else:
                    # for normal building, there is no fragment ion
                    configs[str(record[1])].append(path)

            # delete keys for which there are no valid edge sets
            if len(configs[str(record[1])]) == 0:
                del configs[str(record[1])]

        return configs

    def paths(self, tree, cur=()):
        """
        Parses a tree structure within a dictionary, representing a set of non-isomorphic graphs
        to be used to connect substructures together to generate molecules.

        :param tree: A dictionary containing a set of non-isomorphic graphs for a particular connectivity configuration.

        :param cur: Tuple for results to be appended to.

        :return: For each graph contained within *tree*, generates a tuple of bonds to be formed between substructures.
        """

        if tree == {}:
            yield cur

        else:
            for n, s in tree.items():
                for path in self.paths(s, cur + (n,)):
                    yield path

    def select_substructures(self, l_atoms, table_name):
        """
        Selects specific substructures from a substructure table based on elemental composition. Used to obtain
        sets substructures to be connected.

        :param l_atoms: A list of lists containing the elemental compositions of the queries, in the format
            [C, H, N, O, P, S].

        :param table_name: Name of the substructure table to be queried.

        :return: A list of lists containing the libs of the substructures obtained by the query; the lib is a
            dictionary containing details about the substructure, as generated by
            :py:meth:`metaboblend.databases.get_substructure`.
        """

        subsets = []
        for i in range(len(l_atoms)):

            self.cursor.execute("""SELECT DISTINCT 
                                   smiles,
                                   mol, 
                                   bond_types, 
                                   valence_atoms,  
                                   valence, 
                                   atoms_available,
                                   dummies 
                                       FROM {}
                                       WHERE C = {} 
                                       AND H = {} 
                                       AND N = {} 
                                       AND O = {}
                                       AND P = {}
                                       AND S = {}
                                """.format(table_name, l_atoms[i][0], l_atoms[i][1], l_atoms[i][2],
                                           l_atoms[i][3], l_atoms[i][4], l_atoms[i][5]))
            records = self.cursor.fetchall()
            if len(records) == 0:
                return []

            ss = []
            for record in records:
                ss.append({
                    "smiles": record[0],
                    "mol": Chem.Mol(record[1]),
                    "bond_types": eval(record[2]),
                    "degree_atoms": eval(record[3]),
                    "valence": record[4],
                    "atoms_available": record[5],
                    "dummies": eval(record[6])
                })

            subsets.append(ss)

        return subsets

    def create_compound_database(self):
        """Generates a substructure database, removing previously existing tables if they are present."""

        self.cursor.execute("DROP TABLE IF EXISTS compounds")
        self.cursor.execute("DROP TABLE IF EXISTS substructures")
        self.cursor.execute("DROP TABLE IF EXISTS hmdbid_substructures")
        self.cursor.execute("DROP TABLE IF EXISTS unique_hmdbid")
        self.cursor.execute("DROP TABLE IF EXISTS filtered_hmdbid_substructures")
        self.cursor.execute("DROP TABLE IF EXISTS subset_substructures")

        self.cursor.execute("""CREATE TABLE compounds (
                                   hmdbid TEXT PRIMARY KEY,
                                   exact_mass INTEGER,
                                   formula TEXT,
                                   C INTEGER,
                                   H INTEGER,
                                   N INTEGER,
                                   O INTEGER,
                                   P INTEGER,
                                   S INTEGER,
                                   smiles TEXT)""")

        self.cursor.execute("""CREATE TABLE substructures (
                                   smiles TEXT PRIMARY KEY, 
                                   heavy_atoms INTEGER,
                                   length INTEGER,
                                   exact_mass__1 INTEGER,
                                   exact_mass__0_0001 REAL,
                                   exact_mass REAL,
                                   C INTEGER,
                                   H INTEGER,
                                   N INTEGER,
                                   O INTEGER,
                                   P INTEGER,
                                   S INTEGER,
                                   valence INTEGER,
                                   valence_atoms TEXT,
                                   atoms_available INTEGER,
                                   bond_types TEXT,
                                   dummies TEXT,
                                   mol BLOB)""")

        self.cursor.execute("""CREATE TABLE hmdbid_substructures (
                                   hmdbid TEXT,
                                   smiles,
                                   PRIMARY KEY (hmdbid, smiles))""")

    def create_indexes(self, table="substructures", selection="all"):
        """Creates indexes for the `substructures` table for use by the build method."""

        self.cursor.execute("DROP INDEX IF EXISTS mass__1")
        self.cursor.execute("DROP INDEX IF EXISTS mass__0_0001")
        self.cursor.execute("DROP INDEX IF EXISTS atoms")

        if selection != "gen_subs_table":
            self.cursor.execute("DROP INDEX IF EXISTS heavy_atoms__valence__atoms_available__exact_mass__1")
            self.cursor.execute("DROP INDEX IF EXISTS smiles__heavy_atoms__valence__atoms_available__exact_mass__1")

        self.cursor.execute("""CREATE INDEX mass__1
                               ON %s (exact_mass__1)""" % table)
        self.cursor.execute("""CREATE INDEX mass__0_0001
                               ON %s (exact_mass__0_0001)""" % table)
        self.cursor.execute("""CREATE INDEX atoms 
                               ON %s (C, H, N, O, P, S);""" % table)

        if selection != "gen_subs_table":
            self.cursor.execute("""CREATE INDEX heavy_atoms__valence__atoms_available__exact_mass__1
                                   ON %s (heavy_atoms, atoms_available, valence, exact_mass__1);""" % table)
            self.cursor.execute("""CREATE INDEX smiles__heavy_atoms__valence__atoms_available__exact_mass__1
                                       ON %s (smiles, heavy_atoms, atoms_available, valence, exact_mass__1);""" % table)

    def add_small_substructures(self, table_name):
        """
        Adds a curated list of small substructures (single atom). Similar to update_substructure_database but deals
        with substructures without reference molecules.

        :param table_name: Which table to append the substructure entries to.
        """

        # list of potential relevant small substructures (smarts)
        small_smis = ["CCC", "C=C-C", "C=C=C", "ccc", "ccC", "C=N", "CN", "cnc", "CO", "C=O"]

        for smi in small_smis:
            mol = Chem.MolFromSmarts(smi)

            lib = get_substructure(mol, None)

            smiles_rdkit = Chem.MolToSmiles(lib["mol"])

            exact_mass = calculate_exact_mass(lib["mol"])
            els = get_elements(lib["mol"])

            sub_smi_dict = {'smiles': smiles_rdkit,
                            'exact_mass': exact_mass,
                            'length': sum([els[atom] for atom in els if atom != "*"]),
                            "valence": lib["valence"],
                            "valence_atoms": str(lib["degree_atoms"]),
                            "atoms_available": lib["atoms_available"],
                            "mol": lib["mol"].ToBinary(),
                            "bond_types": str(lib["bond_types"]),
                            "dummies": str(lib["dummies"])}

            sub_smi_dict["exact_mass__1"] = round(sub_smi_dict["exact_mass"], 0)
            sub_smi_dict["exact_mass__0_0001"] = round(sub_smi_dict["exact_mass"], 4)

            sub_smi_dict.update(els)
            sub_smi_dict["heavy_atoms"] = sum([els[atom] for atom in els if atom != "H" and atom != "*"])

            self.cursor.execute("""INSERT OR IGNORE INTO %s (
                                       smiles,
                                       heavy_atoms,
                                       length, 
                                       exact_mass__1, 
                                       exact_mass__0_0001, 
                                       exact_mass, 
                                       C, 
                                       H, 
                                       N, 
                                       O, 
                                       P, 
                                       S, 
                                       valence, 
                                       valence_atoms, 
                                       atoms_available, 
                                       bond_types,
                                       dummies,
                                       mol
                                   ) values (
                                       :smiles,
                                       :heavy_atoms,
                                       :length,
                                       :exact_mass__1,
                                       :exact_mass__0_0001,
                                       :exact_mass,
                                       :C,
                                       :H,
                                       :N,
                                       :O,
                                       :P,
                                       :S,
                                       :valence,
                                       :valence_atoms,
                                       :atoms_available,
                                       :bond_types,
                                       :dummies,
                                       :mol)""" % table_name, sub_smi_dict)

        self.conn.commit()

    def close(self, clean=True):
        if clean:
            self.cursor.execute("DROP TABLE IF EXISTS unique_hmdbid")
            self.cursor.execute("DROP TABLE IF EXISTS filtered_hmdbid_substructures")
            self.cursor.execute("DROP TABLE IF EXISTS subset_substructures")

        self.conn.close()


def get_substructure(mol, idxs_edges_subgraph):
    """
    Generates information for the substructure database from a reference molecule and the bond IDs of a substructure.

    :param mol: An :py:meth:`rdkit.Chem.Mol` object containing a reference molecule that has been fragmented.

    :param idxs_edges_subgraph: List of atom indices within the reference molecule that make up the substructure.

    :return: A list of lists containing the libs of the substructures obtained by the query; the lib is a
        dictionary containing details about the substructure, in the format:

        * "**smiles**": Substructure smiles string

        * "**mol**": Substructure :py:meth:`rdkit.Chem.Mol`

        * "**bond_types**": The type of bonds to be formed by dummy atoms - see
            :py:meth:`metaboblend.build_structures.add_bonds` and :py:meth:`Chem.rdchem.BondType`. Is a dictionary
            whose keys are atom indices and values are bond types, as follows:

            * **1.0** Single
            * **1.5** Aromatic
            * **2.0** Double

        * "**degree_atoms**": A dictionary containing indices of the atoms connected to dummy atoms that can form bonds
            during structure generation as keys, and the number of bonds they can form as values.

        * "**valence**": The total number of bonds that can be formed by the substructure
            (the product of `degree_atoms` and `atoms_available`).

        * "**atoms_available**": The total number of degree atoms.

        * "**dummies**": List of the indices of atoms that may be removed to form bonds during structure generation,
            represented by `*`.
    """

    # convert list of bond indices to list of atom indices
    if idxs_edges_subgraph is None:  # small substructure addition
        atom_idxs_subgraph = [1]
    else:
        atom_idxs_subgraph = []
        for bIdx in idxs_edges_subgraph:
            b = mol.GetBondWithIdx(bIdx)
            a1 = b.GetBeginAtomIdx()
            a2 = b.GetEndAtomIdx()

            if a1 not in atom_idxs_subgraph:
                atom_idxs_subgraph.append(a1)
            if a2 not in atom_idxs_subgraph:
                atom_idxs_subgraph.append(a2)

    # identify atoms which will become dummy elements in the final substructure
    atoms_to_dummy = []
    for idx in atom_idxs_subgraph:
        for atom in mol.GetAtomWithIdx(idx).GetNeighbors():
            if atom.GetIdx() not in atom_idxs_subgraph:
                atoms_to_dummy.append(atom.GetIdx())

    mol_edit = Chem.EditableMol(mol)
    degree_atoms = {}

    for atom in reversed(mol.GetAtoms()):

        if atom.GetIdx() in atoms_to_dummy:
            mol_edit.ReplaceAtom(atom.GetIdx(), Chem.Atom("*"))

    mol = mol_edit.GetMol()
    mol_edit = Chem.EditableMol(mol)

    for atom in reversed(mol.GetAtoms()):
        if atom.GetIdx() not in atom_idxs_subgraph and atom.GetSymbol() != "*":
            mol_edit.RemoveAtom(atom.GetIdx())

    mol_out = mol_edit.GetMol()

    dummies = [atom.GetIdx() for atom in mol_out.GetAtoms() if atom.GetSymbol() == "*"]

    # get bond degrees
    for atom in mol_out.GetAtoms():

        if atom.GetIdx() in dummies:

            for atom_n in atom.GetNeighbors():

                if atom_n.GetSymbol() == "*":
                    continue  # do not count dummies for valence calculations
                elif atom_n.GetIdx() not in degree_atoms:
                    degree_atoms[atom_n.GetIdx()] = 1
                else:
                    degree_atoms[atom_n.GetIdx()] += 1

    # returns the type of the bond as a double (i.e. 1.0 for SINGLE, 1.5 for AROMATIC, 2.0 for DOUBLE)
    bond_types = {}

    for b in mol_out.GetBonds():

        # use bond types to dummy atoms to inform future structure building from compatible substructures
        if mol_out.GetAtomWithIdx(b.GetBeginAtomIdx()).GetSymbol() == "*":
            if b.GetEndAtomIdx() not in bond_types:
                bond_types[b.GetEndAtomIdx()] = [b.GetBondTypeAsDouble()]
            else:
                bond_types[b.GetEndAtomIdx()].append(b.GetBondTypeAsDouble())

        elif mol_out.GetAtomWithIdx(b.GetEndAtomIdx()).GetSymbol() == "*":
            if b.GetBeginAtomIdx() not in bond_types:
                bond_types[b.GetBeginAtomIdx()] = [b.GetBondTypeAsDouble()]
            else:
                bond_types[b.GetBeginAtomIdx()].append(b.GetBondTypeAsDouble())

    try:
        mol_out.UpdatePropertyCache()  # alternative to Chem.SanitizeMol that updates valence information
    except:
        return

    return {"smiles": Chem.MolToSmiles(mol_out),  # REORDERED ATOM INDEXES
            "mol": mol_out,
            "bond_types": bond_types,
            "degree_atoms": degree_atoms,
            "valence": sum(degree_atoms.values()),
            "atoms_available": len(degree_atoms.keys()),
            "dummies": dummies}


def get_elements(mol, elements=None):
    """
    Gets the elemental composition of a molecule.

    :param mol: An :py:meth:`rdkit.Chem.Mol` object containing the molecule of interest.

    :param elements: A dictionary whose keys are strings representing an element and values are 0. Unspecified,
        defaults to `{"C": 0, "H": 0, "N": 0, "O": 0, "P": 0, "S": 0, "*": 0}`.

    :return: The dictionary specified by **elements**, with the number of atoms corresponding to each element as keys.
    """

    if not elements:
        elements = {"C": 0, "H": 0, "N": 0, "O": 0, "P": 0, "S": 0, "*": 0}

    mol = Chem.AddHs(mol)
    for atom in mol.GetAtoms():
        elements[atom.GetSymbol()] += 1

    return elements


def calculate_exact_mass(mol, exact_mass_elements=None):
    """
    Gets the exact mass of a molecule.

    :param mol: An :py:meth:`rdkit.Chem.Mol` object containing the molecule of interest.

    :param exact_mass_elements: A dictionary whose keys are strings representing an element and values are the exact masses of
        each element. Unspecified, defaults to :
        `{"C": 12.0, "H": 1.007825, "N": 14.003074, "O": 15.994915, "P": 30.973763, "S": 31.972072, "*": -1.007825}`

    :return: The exact mass of the molecule.
    """

    if not exact_mass_elements:
        exact_mass_elements = {"C": 12.0, "H": 1.007825, "N": 14.003074, "O": 15.994915, "P": 30.973763, "S": 31.972072,
                               "*": -1.007825}

    exact_mass = 0.0
    mol = Chem.AddHs(mol)
    for atom in mol.GetAtoms():
        atom_symbol = atom.GetSymbol()
        if atom_symbol != "*":
            exact_mass += exact_mass_elements[atom_symbol]

    return exact_mass


def filter_records(records):
    """
    Filters records generated by :py:meth:`parse_xml` to ensure they are compatible with the MetaboBlend workflow.

    :param records: A dictionary containing information about the molecule, as generated by :py:meth:`parse_xml`.

    :return: Generates a dictionary containing key information extracted from the record.

        * "**HMDB_ID**": The HMDB ID of the molecule, in the format "HMDBXXXXXXX".

        * "**formula**": The elemental composition of the molecule, in the format "CXHXNXOXPXSX".

        * "**exact_mass**": The exact mass of the molecule, rounded to 6d.p.

        * "**smiles**": A string containing the smiles representation of the molecule.

        * "**smiles_rdkit**": A string containing the smiles representation of the molecule.

        * "**smiles_rdkit_kek**": A string containing the kekule smiles representation of the molecule.

        * "C", "H", "N", "O", "P", "S": Integers referring to elemental composition.

        * "**mol**": An :py:meth:`rdkit.Chem.Mol` object containing the molecule.
    """

    for record in records:

        if "smiles" in record:
            mol = Chem.MolFromSmiles(record["smiles"])

            try:
                Chem.SanitizeMol(mol)  # confirm mols are legitimate
            except:
                continue

            if mol is None:
                continue

            if mol.GetNumHeavyAtoms() < 4:
                continue  # structures with <3 heavy atoms will not produce useful substructures

            atom_check = [True for atom in mol.GetAtoms() if atom.GetSymbol() not in ["C", "H", "N", "O", "P", "S"]]
            if len(atom_check) > 0:
                continue

            smiles_rdkit = Chem.MolToSmiles(mol)
            smiles_rdkit_kek = Chem.MolToSmiles(mol, kekuleSmiles=True)

            if "+" in smiles_rdkit_kek or "-" in smiles_rdkit_kek or "+" in smiles_rdkit or "-" in smiles_rdkit:
                continue  # only neutral molecules are compatible

            els = get_elements(mol)
            exact_mass = calculate_exact_mass(mol)

            record_dict = {'HMDB_ID': record['accession'],
                           'formula': record["chemical_formula"],
                           'exact_mass': round(exact_mass, 6),
                           'smiles': record['smiles'],
                           'smiles_rdkit': smiles_rdkit,
                           'smiles_rdkit_kek': smiles_rdkit_kek,
                           'C': els['C'],
                           'H': els['H'],
                           'N': els['N'],
                           'O': els['O'],
                           'P': els['P'],
                           'S': els['S'],
                           'mol': mol}

            yield record_dict


def get_substructure_bond_idx(prb_mol, ref_mol):
    """
    Takes a substructure and the original molecule from which it was generated and matches the substructure to its
    original position in the reference molecule. Will only find a single solution, whilst multiple are possible.

    :param prb_mol: Substructure (without dummy atoms) as an :py:meth:`rdkit.Chem.Mol` object.

    :param ref_mol: Original molecule for the substructure to be matched in as an :py:meth:`rdkit.Chem.Mol` object.

    :returns: None if there is no match for the **prb_mol** in the **ref_mol**, else returns a tuple of the bond indices
        in the **ref_mol* matched by the **prb_mol**.
    """

    if ref_mol.HasSubstructMatch(prb_mol):
        atom_idx = ref_mol.GetSubstructMatch(prb_mol)
    else:
        return None

    bond_idx = ()
    for atom in ref_mol.GetAtoms():
        if atom.GetIdx() in atom_idx:
            for bond in atom.GetBonds():
                if bond.GetBeginAtomIdx() in atom_idx and bond.GetEndAtomIdx() in atom_idx:
                    if bond.GetIdx() not in bond_idx:
                        bond_idx = (*bond_idx, bond.GetIdx())

    return bond_idx


def subset_sgs_sizes(sgs, n_min, n_max):
    """
    Some substructure generation methods require that their results be filtered to ensure that they are within the
    databvase substructure size requirements.

    :param sgs: A list of lists containing the indices of edges in the original molecule that represent its
        substructures.

    :param n_min: The minimum number of bonds (edges) for a valid substructure.

    :param n_max: The maximum number of bonds (edges) for a valid substructure.

    :return: The original sgs list, with those substructures that are not within n_min and n_max (inclusive) removed.
    """

    sgs_new = []

    for i, edge_idxs in enumerate(sgs):
        edge_idxs_new = []

        for j, bonds in enumerate(edge_idxs):
            if n_min <= len(bonds) <= n_max:
                edge_idxs_new.append(bonds)

        if len(edge_idxs_new) > 0:
            sgs_new.append(edge_idxs_new)

    return sgs_new


def get_sgs(record_dict, n_min, n_max, method="exhaustive"):
    """
    Generates substructures based on an original molecule, which are described by sets of its bond indices.

    :param record_dict: A dictionary of key information about the original molecule, as generated by
        :py:meth:`metaboblend.databases.filter_records`. Includes HMDBID, smiles and the related
        :py:meth:`rdkit.Chem.Mol` representation.

    :param n_min: The minimum number of bonds (edges) for a valid substructure.

    :param n_max: The maximum number of bonds (edges) for a valid substructure.

    :param method: The method by which to fragment molecules. Substructures must have an exact substructure match in
        the original molecule in order to be considered valid.

        * **exhaustive** The default method for substructure generation. Generates all substructures for a molecule
            within the size range. See :py:meth:`rdkit.Chem.FindAllSubgraphsOfLengthMToN`.

        * **RECAP**  Generates substructures using the retrosynthetic combinatorial analysis procedure; fragments are
            identified that are likely to be useful for drug synthesis. See :py:meth:`rdkit.Chem.RECAP`.

        * **BRICS** Generates substructures by breaking retrosynthetically interesting chemical substructures; fragments
            are identified that are likely to be useful for drug synthesis.. See :py:meth:`rdkit.Chem.BRICS`.

    :return: A list of lists of bond indices referring to substructures of the original molecule.
    """

    if method == "exhaustive":

        return Chem.FindAllSubgraphsOfLengthMToN(record_dict["mol"], n_min, n_max)

    elif method == "RECAP":

        hierarchy = Recap.RecapDecompose(record_dict["mol"])
        sgs = []
        for substructure in hierarchy.GetAllChildren().values():

            substructure = Chem.DeleteSubstructs(substructure.mol, Chem.MolFromSmarts('[#0]'))
            edge_idxs = get_substructure_bond_idx(substructure, record_dict["mol"])
            if edge_idxs is not None:
                sgs.append(edge_idxs)

        return subset_sgs_sizes([sgs], n_min, n_max)

    elif method == "BRICS":

        substructures = BRICS.BRICSDecompose(record_dict["mol"])
        sgs = []
        for substructure in substructures:
            substructure = Chem.DeleteSubstructs(Chem.MolFromSmiles(substructure), Chem.MolFromSmarts('[#0]'))
            edge_idxs = get_substructure_bond_idx(substructure, record_dict["mol"])

            if edge_idxs is not None:
                sgs.append(edge_idxs)

        return subset_sgs_sizes([sgs], n_min, n_max)


def create_substructure_database(hmdb_paths, path_substructure_db, n_min, n_max, method="exhaustive",
                                 max_atoms_available=None, max_valence=None, substructures_only=False):
    """
    Creates a substructure database before using update_substructure_database to add entries. For parameter explanations, see
    update_substructure_database. Also creates indexes on the final database.

    Add entries to the substructure database by fragmenting a set of molecules. Combinations of substructures in this
    database are used to build new molecules.

    :param max_atoms_available: Maximum number of atoms that may be used for bonding by valid substructures.

    :param max_valence: Maximum valence of valid substructures.

    :param hmdb_paths: A list of paths for the HMDB XML record(s) detailing molecules to be fragmented.

    :param path_substructure_db: The path of the SQLite 3 substructure database to be updated.

    :param n_min: The minimum number of bonds (edges) for a valid substructure.

    :param n_max: The maximum number of bonds (edges) for a valid substructure.

    :param records: Records of molecules to be fragmented. Must be a list containing dictionaries containing key
        information about the molecules, as generated by :py:meth:`metaboblend.databases.parse_xml`; if records
        is not supplied, the records will be obtained from the XML at `hmdb_path`.

    :param method: The method by which to fragment molecules. Substructures must have an exact substructure match in
        the original molecule in order to be considered valid.

        * **exhaustive** The default method for substructure generation. Generates all substructures for a molecule
            within the size range. See :py:meth:`rdkit.Chem.FindAllSubgraphsOfLengthMToN`.

        * **RECAP**  Generates substructures using the retrosynthetic combinatorial analysis procedure; fragments are
            identified that are likely to be useful for drug synthesis. See :py:meth:`rdkit.Chem.RECAP`.

        * **BRICS** Generates substructures by breaking retrosynthetically interesting chemical substructures; fragments
            are identified that are likely to be useful for drug synthesis.. See :py:meth:`rdkit.Chem.BRICS`.

    :param substructures_only: Whether to generate all tables or only the substructures table.
    """

    db = SubstructureDb(path_substructure_db)
    db.create_compound_database()
    db.close()

    for hmdb_path in hmdb_paths:
        update_substructure_database(hmdb_path, path_substructure_db, n_min, n_max, method=method,
                                     max_atoms_available=max_atoms_available, max_valence=max_valence,
                                     substructures_only=substructures_only)

    db = SubstructureDb(path_substructure_db)
    db.create_indexes()
    db.close()


def update_substructure_database(hmdb_path, path_substructure_db, n_min, n_max, records=None, method="exhaustive",
                                 max_atoms_available=None, max_valence=None, substructures_only=False):
    """
    Add entries to the substructure database by fragmenting a set of molecules. Combinations of substructures in this
    database are used to build new molecules.

    :param max_atoms_available: Maximum number of atoms that may be used for bonding by valid substructures.

    :param max_valence: Maximum valence of valid substructures.

    :param hmdb_path: The path of the HMDB XML record(s) detailing molecules to be fragmented. Will be overriden by
        `records` if provided.

    :param path_substructure_db: The path of the SQLite 3 substructure database to be updated.

    :param n_min: The minimum number of bonds (edges) for a valid substructure.

    :param n_max: The maximum number of bonds (edges) for a valid substructure.

    :param records: Records of molecules to be fragmented. Must be a list containing dictionaries containing key
        information about the molecules, as generated by :py:meth:`metaboblend.databases.parse_xml`; if records
        is not supplied, the records will be obtained from the XML at `hmdb_path`.

    :param method: The method by which to fragment molecules. Substructures must have an exact substructure match in
        the original molecule in order to be considered valid.

        * **exhaustive** The default method for substructure generation. Generates all substructures for a molecule
            within the size range. See :py:meth:`rdkit.Chem.FindAllSubgraphsOfLengthMToN`.

        * **RECAP**  Generates substructures using the retrosynthetic combinatorial analysis procedure; fragments are
            identified that are likely to be useful for drug synthesis. See :py:meth:`rdkit.Chem.RECAP`.

        * **BRICS** Generates substructures by breaking retrosynthetically interesting chemical substructures; fragments
            are identified that are likely to be useful for drug synthesis.. See :py:meth:`rdkit.Chem.BRICS`.

    :param substructures_only: Whether to generate all tables or only the substructures table.
    """

    conn = sqlite3.connect(path_substructure_db)
    cursor = conn.cursor()

    if records is None:
        records = parse_xml(hmdb_path, reformat=False)


    for record_dict in filter_records(records):
        if not substructures_only:
            cursor.execute("""INSERT OR IGNORE INTO compounds (
                                   hmdbid, 
                                   exact_mass, 
                                   formula, 
                                   C, H, N, O, P, S, 
                                   smiles)
                               VALUES (
                                   :HMDB_ID, 
                                   :exact_mass,
                                   :formula, 
                                   :C, :H, :N, :O, :P, :S, 
                                   :smiles)""", record_dict)

        # Returns a tuple of 2-tuples with bond IDs
        for sgs in get_sgs(record_dict, n_min, n_max, method=method):
            for edge_idxs in sgs:
                lib = get_substructure(record_dict["mol"], edge_idxs)  # convert bond IDs to substructure mol

                if lib is None:
                    continue

                if max_atoms_available is not None:
                    if lib["atoms_available"] > max_atoms_available:
                        continue

                if max_valence is not None:
                    if lib["valence"] > max_valence:
                        continue

                smiles_rdkit = Chem.MolToSmiles(lib["mol"])  # canonical rdkit smiles

                exact_mass = calculate_exact_mass(lib["mol"])
                els = get_elements(lib["mol"])

                sub_smi_dict = {'smiles': smiles_rdkit,
                                'exact_mass': exact_mass,
                                'length': sum([els[atom] for atom in els if atom != "*"]),
                                "valence": lib["valence"],
                                "valence_atoms": str(lib["degree_atoms"]),
                                "atoms_available": lib["atoms_available"],
                                "mol": lib["mol"].ToBinary(),
                                "bond_types": str(lib["bond_types"]),
                                "dummies": str(lib["dummies"])}

                sub_smi_dict["exact_mass__1"] = round(sub_smi_dict["exact_mass"], 0)
                sub_smi_dict["exact_mass__0_0001"] = round(sub_smi_dict["exact_mass"], 4)

                sub_smi_dict.update(els)
                sub_smi_dict["heavy_atoms"] = sum([els[atom] for atom in els if atom != "H" and atom != "*"])

                cursor.execute("""INSERT OR IGNORE INTO substructures (
                                      smiles, 
                                      heavy_atoms, 
                                      length, 
                                      exact_mass__1, 
                                      exact_mass__0_0001, 
                                      exact_mass, 
                                      C, 
                                      H, 
                                      N, 
                                      O, 
                                      P, 
                                      S, 
                                      valence, 
                                      valence_atoms, 
                                      atoms_available, 
                                      bond_types,
                                      dummies,
                                      mol)
                                  values (
                                      :smiles,
                                      :heavy_atoms,
                                      :length,
                                      :exact_mass__1,
                                      :exact_mass__0_0001,
                                      :exact_mass,
                                      :C,
                                      :H,
                                      :N,
                                      :O,
                                      :P,
                                      :S,
                                      :valence,
                                      :valence_atoms,
                                      :atoms_available,
                                      :bond_types,
                                      :dummies,
                                      :mol)""", sub_smi_dict)

                if not substructures_only:
                    cursor.execute("""INSERT OR IGNORE INTO hmdbid_substructures (
                                          hmdbid, 
                                          smiles) 
                                      VALUES ("%s", "%s")""" % (record_dict['HMDB_ID'], smiles_rdkit))

    conn.commit()
    conn.close()


def create_isomorphism_database(db_out, max_n_substructures, max_atoms_available, path_ri=None):
    """
    Generates a connectivity database containing sets of possible combinations of substructures; also generates PKL
    files containing the graphs required for building molecules from substructures. The connectivity database is
    required to build molecules from substructures.

    :param db_out: The path of the SQLite 3 database to be generated.

    :param max_n_substructures: The maximal number of substructures (vertices). At least two substructures must be
        available for bonding for a graph to be created.

    :param max_atoms_available: The maximal number of atoms available (maximal number of edges per vertice) in each
        substructure for bonding. At least one atom must be available for bonding for a graph to be created.

    :param path_ri: The path of RI, a tool for verifying subgraph isomorphism.
    """
    
    conn = sqlite3.connect(db_out)
    cursor = conn.cursor()

    cursor.execute("""DROP TABLE IF EXISTS subgraphs""")
    cursor.execute("""CREATE TABLE subgraphs (
                          id_pkl INTEGER,
                          n_graphs INTEGER,
                          graph6 TEXT,
                          k INTEGER,
                          k_partite TEXT,
                          k_valences TEXT,
                          nodes_valences TEXT,
                          n_nodes INTEGER,
                          n_edges INTEGER,
                          root BLOB,
                          PRIMARY KEY (graph6, k_partite, nodes_valences)
                   );""")
    conn.commit()

    id_pkl = 0

    for g, p in calculate_complete_multipartite_graphs(max_atoms_available, max_n_substructures):

        # get complete set of non-isomorphic graphs, using geng, from a distinct multipartite graph as input
        proc = subprocess.Popen(["geng", str(len(g.nodes)), "-d1", "-D2", "-q"], stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)  # max valence for single atom of 2
        geng_out, err = proc.communicate()

        proc.stdout.close()
        proc.stderr.close()

        # pipe geng output to RI to generate mappings (complete set of non-isomorphic configurations)
        for i, line_geng in enumerate(geng_out.split()):

            s_g = nx.read_graph6(io.BytesIO(line_geng))

            k_gfu = tempfile.NamedTemporaryFile(mode="w", delete=False)
            k_gfu.write(graph_to_ri(g, "k_graph"))
            k_gfu.seek(0)

            s_gfu = tempfile.NamedTemporaryFile(mode="w", delete=False)
            s_gfu.write(graph_to_ri(s_g, "subgraph"))
            s_gfu.seek(0)

            proc = subprocess.Popen([path_ri, "mono", "geu", k_gfu.name, s_gfu.name], stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            ri_out, err = proc.communicate()

            k_gfu.close()
            s_gfu.close()

            mappings = []
            subgraphs = {}

            for line in ri_out.decode("utf-8").splitlines():
                if line[0] == "{":
                    mappings.append(eval(line))

            if len(mappings) > 0:
                gi = graph_info(p, s_g, mappings, )  # convert mappings to valence/connectivity specifications

                for vn in gi:
                    if vn not in subgraphs:
                        subgraphs[vn] = gi[vn]

                    else:
                        for es in gi[vn]:
                            if es not in subgraphs[vn]:
                                subgraphs[vn].append(es)

            if len(subgraphs) > 0:
                for vn in subgraphs:  # for each valence configuration
                    subgraphs[vn] = sort_subgraphs(subgraphs[vn])  # sort to remove duplicate configurations
                    root = {}  # graph to be pickled

                    for fr in subgraphs[vn]:
                        parent = root
                        for e in fr:
                            parent = parent.setdefault(e, {})

                    vt = tuple([sum(v) for v in eval(vn)])

                    id_pkl += 1
                    cursor.execute("""INSERT INTO subgraphs (
                                          id_pkl, 
                                          n_graphs, 
                                          graph6,
                                          k,
                                          k_partite,
                                          k_valences,
                                          nodes_valences,
                                          n_nodes, n_edges,
                                          root)
                                      VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", (
                                          id_pkl,
                                          len(subgraphs[vn]),
                                          line_geng,
                                          len(p),
                                          str(p),
                                          str(vt),
                                          str(vn),
                                          s_g.number_of_nodes(),
                                          s_g.number_of_edges(),
                                          pickle.dumps(root)))

            conn.commit()
    conn.close()
