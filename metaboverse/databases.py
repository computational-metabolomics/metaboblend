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

import io
import os
import sys
import subprocess
import pickle
import sqlite3
import tempfile
from collections import OrderedDict
import xml.etree.ElementTree as etree
import networkx as nx
from rdkit import Chem
from rdkit.Chem import Recap
from rdkit.Chem import BRICS
from .auxiliary import calculate_complete_multipartite_graphs, graph_to_ri, graph_info, sort_subgraphs, draw_subgraph

sqlite3.register_converter("PICKLE", pickle.loads)


def reformat_xml(source, encoding="utf8"):
    """
    Reformats HMDB xml files to be compatible with :py:meth:`metaboverse.databases.parse_xml`; some such files do not
    contain a `<hmdb xmlns="http://www.hmdb.ca">` header.

    :param source: File to be reformatted.

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

    :param reformat: Whether to apply :py:meth:`metaboverse.databases.reformat_xml` to the XML file. Is required for
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

                for event, elem in etree.iterparse(inp, events=("start", "end")):
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


class ConnectivityDb:
    """
    Object containing a reference to the connectivity database.

    :param db: Path to the connectivity database.
    """

    def __init__(self, db):
        """Constructor method"""

        self.db = db


class SubstructureDb:
    """
    Methods for interacting with the SQLITE3 substructure and connectivity databases. Provides a connection to the
    substructure database and, if provided, the connectivity database.

    :param db: Path to the substructure database.
    :param db2: Path to the connectivity database.
    :ivar path_pkls: Path to the directory containing connectivity graph PKLs
    :param conn: A :py:meth:`sqlite3.connection` to the substructure database; the connectivity database will be attached
        as 'graphs'.
    :param cursor: A :py:meth:`sqlite3.connection.cursor` for the substructure database; the connectivity database will
        be attached as "graphs".
    """

    def __init__(self, db, path_pkls="", db2=None):
        """Constructor method"""

        self.db = db
        self.db2 = db2
        self.path_pkls = path_pkls

        self.conn = sqlite3.connect(self.db)
        self.cursor = self.conn.cursor()

        if self.db2 is not None:
            self.cursor.execute("ATTACH DATABASE '%s' as 'graphs';" % self.db2)

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

        self.cursor.execute("""SELECT DISTINCT hmdbid, exact_mass, formula, C, H, N, O, P, S, smiles,
                               smiles_rdkit, smiles_rdkit_kek FROM compounds%s""" % sql)

        return self.cursor.fetchall()

    def filter_hmdbid_substructures(self, min_node_weight=2):
        """
        Filters `hmdbid_substructures` table in the substructure database by the frequency of the substructure; allows
        for the removal of non-unique (min_node_weight = 2) or infrequent (min_node_weight > 2) substructures. Generates
        a new table, `unique_hmdbid`, which contains distinct hmdbids and 'filtered_hmdbid_substructures' which contains
        the filtered substructures with their frequency.

        :param min_node_weight: Minimal count of the substructure within 'hmdbid_substructures'.
        """

        self.cursor.execute("DROP TABLE IF EXISTS unique_hmdbid")
        self.cursor.execute("DROP TABLE IF EXISTS filtered_hmdbid_substructures")

        self.cursor.execute("CREATE TABLE unique_hmdbid AS SELECT DISTINCT hmdbid FROM compounds")

        self.cursor.execute("""CREATE TABLE filtered_hmdbid_substructures AS
                                   SELECT smiles_rdkit, COUNT(*) FROM hmdbid_substructures
                                   GROUP BY smiles_rdkit HAVING COUNT(*) >=%s""" % min_node_weight)

    def generate_substructure_network(self, method="default", min_node_weight=2, remove_isolated=False):
        """
        Generate networks to explore the co-occurence of substructures in the substructure database.
        
        :param method: Which method to use for the generation of a substructure co-occurence network.

            * **default** Generate a standard substructure network using the most time-efficient method. Node weights
                represent the frequency of substructures in the substructure database whilst edge weights represent the
                co-occurence of these substructures. See
                :py:meth:`metaboverse.databases.default_substructure_network`.

            * **extended** Alternative method for generating substructure networks that generates an intermediate
                metabolite-substructure linkage network; returns the same network as the **default** method. See
                :py:meth:`metaboverse.databases.extended_substructure_network`.

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
        generated by :py:meth:`metaboverse.databases.default_substructure_network`.

        :param substructure_graph: A :py:meth:`networkx.Graph` containing substructure nodes, weighted by their
            frequencies in the substructure database. Can be passed by
            :py:meth:`metaboverse.databases.generate_substructure_network`.

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
        self.cursor.execute("""SELECT * FROM hmdbid_substructures WHERE smiles_rdkit IN 
                                   (SELECT smiles_rdkit FROM filtered_hmdbid_substructures)""")
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
            :py:meth:`metaboverse.databases.generate_substructure_network`.

        :param unique_hmdb_ids: A list containing distinct HMDB IDs from the substructure database.

        :return: A :py:meth:`networkx.Graph` object containing the specified network.
        """

        # add edges by walking through hmdbid_substructures
        for unique_hmdb_id in unique_hmdb_ids:
            self.cursor.execute("""SELECT * FROM hmdbid_substructures 
                                       WHERE smiles_rdkit IN (SELECT smiles_rdkit FROM filtered_hmdbid_substructures) 
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
            * **0_1** One decimal place
            * **0_01** Two decimal places
            * **0_001** Three decimal places
            * **0_0001** Four decimal places

        :param masses: A list of integers to limit the query by.

        :param table_name: Name of the substructure table to be queried.

        :return: Sorted list of mass values from the substructure database, filtered by the supplied parameters.
        """

        mass_values = []
        filter_mass = ""
        if type(masses) == list:
            if len(masses) > 0:
                filter_mass = " WHERE exact_mass__1 IN ({})".format(",".join(map(str, masses)))

        self.cursor.execute("SELECT DISTINCT exact_mass__{} FROM {}{}".format(accuracy, table_name, filter_mass))

        records = self.cursor.fetchall()
        for record in records:
            mass_values.append(record[0])
        mass_values.sort()

        return mass_values

    def select_ecs(self, exact_mass, table_name, accuracy, ppm=None):
        """
        Select elemental compositions based on an exact mass; allows for inexact mass searches based on error (ppm).

        :param exact_mass: The exact mass used to query the substructure database.

        :param table_name: Name of the substructure table to be queried.

        :param accuracy: To which decimal places of accuracy results are to be limited to.

            * **1** Integer level
            * **0_1** One decimal place
            * **0_01** Two decimal places
            * **0_001** Three decimal places
            * **0_0001** Four decimal places

        :param ppm: The allowable error of the query (in parts per million).

        :return: Returns the elemental compositions associated with the `exact_mass` via
            :py:meth:`sqlite3.connection.cursor.fetchall`; essentially provides a list containing a list for the values
            of each row (C, H, N, O, P, S).
        """

        if ppm is None:
            mass_statement = "= " + str(exact_mass)
        else:
            tolerance = (exact_mass / 1000000) * ppm
            mass_statement = "< {} AND exact_mass__{} > {}".format(exact_mass + tolerance,
                                                                   accuracy,
                                                                   exact_mass - tolerance)

        self.cursor.execute("""SELECT DISTINCT C, H, N, O, P, S
                                   FROM {} 
                                   WHERE exact_mass__{} {}
                            """.format(table_name, accuracy, mass_statement))

        return self.cursor.fetchall()

    def paths(self, tree, cur=()):
        if tree == {}:
            yield cur
        else:
            for n, s in tree.items():
                for path in self.paths(s, cur + (n,)):
                    yield path

    def isomorphism_graphs(self, id_pkl):
        with open(os.path.join(self.path_pkls, "{}.pkl".format(id_pkl)), 'rb') as pickle_file:
            nGcomplete = pickle.load(pickle_file)
        for p in self.paths(nGcomplete):
            yield p

    def k_configs(self):
        """
        Obtains strings detailing the valences for each substructure in a connectivity graph and the ID of the related
        PKL file. Used to match a set of substructures to the correct set of non-isomorphic graphs in the connectivity
        database.

        :return: Dictionary containing the valences as keys and PKL IDs as values.
        """

        self.cursor.execute("""SELECT id_pkl, nodes_valences 
                                   FROM subgraphs""")

        records = self.cursor.fetchall()
        configs = {}

        for record in records:
            configs[str(record[1])] = record[0]

        return configs

    def select_sub_structures(self, l_atoms, table_name):
        """
        Selects specific substructures from a substructure table based on elemental composition. Used to obtain
        sets substructures to be connected.

        :param l_atoms: A list of lists containing the elemental compositions of the queries, in the format
            [C, H, N, O, P, S].

        :param table_name: Name of the substructure table to be queried.

        :return: A list of lists containing the libs of the substructures obtained by the query; the lib is a
            dictionary containing details about the substructure, as generated by
            :py:meth:`metaboverse.databases.get_substructure`.
        """

        subsets = []
        for i in range(len(l_atoms)):

            self.cursor.execute("""SELECT DISTINCT lib 
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
            ss = [pickle.loads(record[0]) for record in records]
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
                                   smiles TEXT,
                                   smiles_rdkit TEXT,
                                   smiles_rdkit_kek TEXT)""")

        self.cursor.execute("""CREATE TABLE substructures (
                                   smiles TEXT PRIMARY KEY, 
                                   heavy_atoms INTEGER,
                                   length INTEGER,
                                   exact_mass__1 INTEGER,
                                   exact_mass__0_1 REAL,
                                   exact_mass__0_01 REAL,
                                   exact_mass__0_001 REAL,
                                   exact_mass__0_0001 REAL,
                                   exact_mass REAL,
                                   count INTEGER,
                                   C INTEGER,
                                   H INTEGER,
                                   N INTEGER,
                                   O INTEGER,
                                   P INTEGER,
                                   S INTEGER,
                                   valence INTEGER,
                                   valence_atoms TEXT,
                                   atoms_available INTEGER,
                                lib PICKLE)""")

        self.cursor.execute("""CREATE TABLE hmdbid_substructures (
                                   hmdbid TEXT,
                                   smiles_rdkit,
                                   PRIMARY KEY (hmdbid, smiles_rdkit))""")

    def create_indexes(self):
        """Creates indexes for the `substructures` table for use by the build method."""

        self.cursor.execute("DROP INDEX IF EXISTS heavy_atoms__Valence__mass__1__idx")
        self.cursor.execute("DROP INDEX IF EXISTS heavy_atoms__Valence__mass__0_1__idx")
        self.cursor.execute("DROP INDEX IF EXISTS heavy_atoms__Valence__mass__0_01__idx")
        self.cursor.execute("DROP INDEX IF EXISTS heavy_atoms__Valence__mass__0_001__idx")
        self.cursor.execute("DROP INDEX IF EXISTS heavy_atoms__Valence__mass__0_0001__idx")
        self.cursor.execute("DROP INDEX IF EXISTS atoms__Valence__idx")

        self.cursor.execute("""CREATE INDEX heavy_atoms__Valence__mass__1__idx 
                               ON substructures (heavy_atoms, valence, valence_atoms, exact_mass__1)""")
        self.cursor.execute("""CREATE INDEX heavy_atoms__Valence__mass__0_1__idx 
                               ON substructures (heavy_atoms, valence, valence_atoms, exact_mass__0_1)""")
        self.cursor.execute("""CREATE INDEX heavy_atoms__Valence__mass__0_01__idx 
                               ON substructures (heavy_atoms, valence, valence_atoms, exact_mass__0_01)""")
        self.cursor.execute("""CREATE INDEX heavy_atoms__Valence__mass__0_001__idx 
                               ON substructures (heavy_atoms, valence, valence_atoms, exact_mass__0_001)""")
        self.cursor.execute("""CREATE INDEX heavy_atoms__Valence__mass__0_0001__idx
                               ON substructures (heavy_atoms, valence, valence_atoms, exact_mass__0_0001)""")
        self.cursor.execute("""CREATE INDEX atoms__Valence__idx 
                               ON substructures (C, H, N, O, P, S, valence, valence_atoms);""")

    def close(self):
        self.cursor.execute("DROP TABLE IF EXISTS unique_hmdbid")
        self.cursor.execute("DROP TABLE IF EXISTS filtered_hmdbid_substructures")
        self.cursor.execute("DROP TABLE IF EXISTS subset_substructures")
        self.conn.close()


def get_substructure(mol, idxs_edges_subgraph, debug=False):
    """
    Generates information for the substructure database from a reference molecule and the bond IDs of a substructure.

    :param mol: An :py:meth:`rdkit.Chem.Mol' object containing a reference molecule that has been fragmented.

    :param idxs_edges_subgraph: List of atom indices within the reference molecule that make up the substructure.

    :return: A list of lists containing the libs of the substructures obtained by the query; the lib is a
        dictionary containing details about the substructure, in the format:

        * "**smiles**": Substructure smiles string

        * "**mol**": Substructure :py:meth:`rdkit.Chem.Mol`

        * "**bond_types**": The type of bonds to be formed by dummy atoms - see
            :py:meth:`metaboverse.build_structures.add_bonds` and :py:meth:`Chem.rdchem.BondType`. Is a dictionary
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

    atom_idxs_subgraph = []
    for bIdx in idxs_edges_subgraph:
        b = mol.GetBondWithIdx(bIdx)
        a1 = b.GetBeginAtomIdx()
        a2 = b.GetEndAtomIdx()

        if a1 not in atom_idxs_subgraph:
            atom_idxs_subgraph.append(a1)
        if a2 not in atom_idxs_subgraph:
            atom_idxs_subgraph.append(a2)

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

    for atom in mol_out.GetAtoms():

        if atom.GetIdx() in dummies:

            for atom_n in atom.GetNeighbors():

                if atom_n.GetSymbol() == "*":
                    continue  # do not count dummies for valence calculations
                elif atom_n.GetIdx() not in degree_atoms:
                    degree_atoms[atom_n.GetIdx()] = 1
                else:
                    degree_atoms[atom_n.GetIdx()] += 1

    # Returns the type of the bond as a double (i.e. 1.0 for SINGLE, 1.5 for AROMATIC, 2.0 for DOUBLE)
    bond_types = {}

    for b in mol_out.GetBonds():
        if debug:
            print(b.GetBondTypeAsDouble())
            print(b.GetBondType())
            print(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), mol_out.GetAtomWithIdx(b.GetBeginAtomIdx()).GetSymbol(),
                  mol_out.GetAtomWithIdx(b.GetEndAtomIdx()).GetSymbol())

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
        atomSymbol = atom.GetSymbol()
        if atomSymbol != "*":
            exact_mass += exact_mass_elements[atomSymbol]

    return exact_mass


def filter_records(records):
    """
    Filters records generated by :py:meth:`parse_xml` to ensure they are compatible with the Metaboverse workflow.

    :param records: A dictionary containing information about the molecule, as generated by :py:meth:`parse_xml`.

    :return: Generates a dictionary containing key information extracted from the record.

        * "**HMDB_ID**": The HMDB ID of the molecule, in the format "HMDBXXXXXXX".

        * "**formula**": The elemental composition of the molecule, in the format "CXHXNXOXPXSX".

        * "**exact_mass**": The exact mass of the molecule, rounded to 6d.p.

        * "**smiles**": A string containing the smiles representation of the molecule.

        * "**smiles_rdkit**": A string containing the smiles representation of the molecule.

        * "**smiles_rdkit_kek**": A string containing the kekule smiles representation of the molecule.

        * "C", "H", "N", "O", "P", "S": Integers referring to elemental composition.

        * "**mol**": An :py:meth:`rdkit.Chem.Mol' object containing the molecule.
    """

    for record in records:

        if "smiles" in record:
            mol = Chem.MolFromSmiles(record["smiles"])
            try:
                Chem.SanitizeMol(mol)
            except:
                continue

            if mol is None:
                continue

            if mol.GetNumHeavyAtoms() < 4:
                continue

            atom_check = [True for atom in mol.GetAtoms() if atom.GetSymbol() not in ["C", "H", "N", "O", "P", "S"]]
            if len(atom_check) > 0:
                continue

            smiles_rdkit = Chem.MolToSmiles(mol)
            smiles_rdkit_kek = Chem.MolToSmiles(mol, kekuleSmiles=True)

            if "+" in smiles_rdkit_kek or "-" in smiles_rdkit_kek or "+" in smiles_rdkit or "-" in smiles_rdkit:
                continue

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


def subset_sgs_max_atoms_available(sgs, n_min, n_max):
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
        :py:meth:`metaboverse.databases.filter_records`. Includes HMDBID, smiles and the related
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

        return subset_sgs_max_atoms_available([sgs], n_min, n_max)

    elif method == "BRICS":

        substructures = BRICS.BRICSDecompose(record_dict["mol"])
        sgs = []
        for substructure in substructures:
            substructure = Chem.DeleteSubstructs(Chem.MolFromSmiles(substructure), Chem.MolFromSmarts('[#0]'))
            edge_idxs = get_substructure_bond_idx(substructure, record_dict["mol"])

            if edge_idxs is not None:
                sgs.append(edge_idxs)

        return subset_sgs_max_atoms_available([sgs], n_min, n_max)


def update_substructure_database(fn_hmdb, fn_db, n_min, n_max, records=None, method="exhaustive"):
    """
    Add entries to the substructure database by fragmenting a set of molecules. Combinations of substructures in this
    database are used to build new molecules.

    :param fn_hmdb: The path of the HMDB XML record(s) detailing molecules to be fragmented. Will be overriden by
        `records` if provided.

    :param fn_db: The path of the SQLite 3 substructure database to be updated.

    :param n_min: The minimum number of bonds (edges) for a valid substructure.

    :param n_max: The maximum number of bonds (edges) for a valid substructure.

    :param records: Records of molecules to be fragmented. Must be a list containing dictionaries containing key
        information about the molecules, as generated by :py:meth:`metaboverse.databases.parse_xml`; if records
        is not supplied, the records will be obtained from the XML at `fn_hmdb`.

    :param method: The method by which to fragment molecules. Substructures must have an exact substructure match in
        the original molecule in order to be considered valid.

        * **exhaustive** The default method for substructure generation. Generates all substructures for a molecule
            within the size range. See :py:meth:`rdkit.Chem.FindAllSubgraphsOfLengthMToN`.

        * **RECAP**  Generates substructures using the retrosynthetic combinatorial analysis procedure; fragments are
            identified that are likely to be useful for drug synthesis. See :py:meth:`rdkit.Chem.RECAP`.

        * **BRICS** Generates substructures by breaking retrosynthetically interesting chemical substructures; fragments
            are identified that are likely to be useful for drug synthesis.. See :py:meth:`rdkit.Chem.BRICS`.
    """

    conn = sqlite3.connect(fn_db)
    cursor = conn.cursor()

    if records is None:
        records = parse_xml(fn_hmdb, reformat=False)

    for record_dict in filter_records(records):

        cursor.execute("""INSERT OR IGNORE INTO compounds (
                               hmdbid, 
                               exact_mass, 
                               formula, 
                               C, H, N, O, P, S, 
                               smiles, 
                               smiles_rdkit, 
                               smiles_rdkit_kek)
                           VALUES (
                               :HMDB_ID, 
                               :exact_mass,
                               :formula, 
                               :C, :H, :N, :O, :P, :S, 
                               :smiles, 
                               :smiles_rdkit, 
                               :smiles_rdkit_kek)""", record_dict)

        # Returns a tuple of 2-tuples with bond IDs
        for sgs in get_sgs(record_dict, n_min, n_max, method=method):
            for edge_idxs in sgs:
                lib = get_substructure(record_dict["mol"], edge_idxs)
                if lib is None:
                    continue

                smiles_rdkit = Chem.MolToSmiles(lib["mol"])

                exact_mass = calculate_exact_mass(lib["mol"])
                els = get_elements(lib["mol"])

                pkl_lib = pickle.dumps(lib)
                sub_smi_dict = {'smiles': smiles_rdkit,
                                'exact_mass': exact_mass,
                                'count': 0,
                                'length': sum([els[atom] for atom in els if atom != "*"]),
                                "valence": lib["valence"],
                                "valence_atoms": str(lib["degree_atoms"]),
                                "atoms_available": lib["atoms_available"],
                                "lib": pkl_lib}

                sub_smi_dict["exact_mass__1"] = round(sub_smi_dict["exact_mass"], 0)
                sub_smi_dict["exact_mass__0_1"] = round(sub_smi_dict["exact_mass"], 1)
                sub_smi_dict["exact_mass__0_01"] = round(sub_smi_dict["exact_mass"], 2)
                sub_smi_dict["exact_mass__0_001"] = round(sub_smi_dict["exact_mass"], 3)
                sub_smi_dict["exact_mass__0_0001"] = round(sub_smi_dict["exact_mass"], 4)

                sub_smi_dict.update(els)
                sub_smi_dict["heavy_atoms"] = sum([els[atom] for atom in els if atom != "H" and atom != "*"])

                cursor.execute("""INSERT OR IGNORE INTO substructures (
                                      smiles, 
                                      heavy_atoms, 
                                      length, 
                                      exact_mass__1, 
                                      exact_mass__0_1, 
                                      exact_mass__0_01, 
                                      exact_mass__0_001, 
                                      exact_mass__0_0001, 
                                      exact_mass, count, 
                                      C, 
                                      H, 
                                      N, 
                                      O, 
                                      P, 
                                      S, 
                                      valence, 
                                      valence_atoms, 
                                      atoms_available, 
                                      lib)
                                  values (
                                      :smiles,
                                      :heavy_atoms,
                                      :length,
                                      :exact_mass__1,
                                      :exact_mass__0_1,
                                      :exact_mass__0_01,
                                      :exact_mass__0_001,
                                      :exact_mass__0_0001,
                                      :exact_mass,
                                      :count,
                                      :C,
                                      :H,
                                      :N,
                                      :O,
                                      :P,
                                      :S,
                                      :valence,
                                      :valence_atoms,
                                      :atoms_available,:lib)""", sub_smi_dict)

                cursor.execute("""INSERT OR IGNORE INTO hmdbid_substructures (
                                      hmdbid, 
                                      smiles_rdkit) 
                                  VALUES ("%s", "%s")""" % (record_dict['HMDB_ID'], smiles_rdkit))
    conn.commit()
    conn.close()


def create_isomorphism_database(db_out, pkls_out, max_n_substructures, max_atoms_available, path_geng=None, path_RI=None, debug=False):
    """
    Generates a connectivity database containing sets of possible combinations of substructures; also generates PKL
    files containing the graphs required for building molecules from substructures. The connectivity database is
    required to build moelcules from substructures.

    :param db_out: The path of the SQLite 3 database to be generated.

    :param pkls_out: The directory in which to dump the PKL files containing teh graphs required for building 
        molecules from substructures.

    :param max_n_substructures: The maximal number of substructures (vertices). At least two substructures must be
        available for bonding for a graph to be created.

    :param max_atoms_available: The maximal number of atoms available (maximal number of edges per vertice) in each
        substructure for bonding. At least one atom must be available for bonding for a graph to be created.

    :param path_geng: The path of geng, a tool for the generation of small non-isomorphic graphs.

    :param path_RI: The path of RI, a tool for verifying subgraph isomorphism.

    :param debug: Debug print statements provide further information on how the function is generating the connectivity
        database.
        * **True** Print debug statements.
        * **False** Hide debug print statements.
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
                          PRIMARY KEY (graph6, k_partite, nodes_valences)
                   );""")
    conn.commit()

    id_pkl = 0

    for G, p in calculate_complete_multipartite_graphs(max_atoms_available, max_n_substructures):

        if debug:
            print([path_geng, str(G.number_of_nodes()), "-d1", "-D2", "-q"])  # max valence for single atom of 2
        proc = subprocess.Popen([path_geng, str(len(G.nodes)), "-d1", "-D2", "-q"], stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        geng_out, err = proc.communicate()

        proc.stdout.close()
        proc.stderr.close()

        for i, line_geng in enumerate(geng_out.split()):
            
            if debug:
                print(line_geng)

            sG = nx.read_graph6(io.BytesIO(line_geng))

            k_gfu = tempfile.NamedTemporaryFile(mode="w", delete=False)
            k_gfu.write(graph_to_ri(G, "k_graph"))
            k_gfu.seek(0)

            s_gfu = tempfile.NamedTemporaryFile(mode="w", delete=False)
            s_gfu.write(graph_to_ri(sG, "subgraph"))
            s_gfu.seek(0)

            proc = subprocess.Popen([path_RI, "mono", "geu", k_gfu.name, s_gfu.name], stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            RI_out, err = proc.communicate()

            k_gfu.close()
            s_gfu.close()

            mappings = []
            subgraphs = {}

            for line in RI_out.decode("utf-8").splitlines():
                if line[0] == "{":
                    mappings.append(eval(line))

            if len(mappings) > 0:
                gi = graph_info(p, sG, mappings, )

                for vn in gi:
                    if vn not in subgraphs:
                        subgraphs[vn] = gi[vn]

                    else:
                        for es in gi[vn]:
                            if es not in subgraphs[vn]:
                                subgraphs[vn].append(es)

            if len(subgraphs) > 0:
                for vn in subgraphs:
                    subgraphs[vn] = sort_subgraphs(subgraphs[vn])
                    root = {}

                    if debug:
                        col_plt, sG_plt = draw_subgraph(subgraphs[vn][0], eval(vn))
                        col_plt.show()

                    for fr in subgraphs[vn]:
                        parent = root
                        for e in fr:
                            parent = parent.setdefault(e, {})

                    vt = tuple([sum(v) for v in eval(vn)])
                    if debug:
                        print("INSERT:", i, line_geng.decode("utf-8"), len(subgraphs[vn]), len(p), str(p), vt, vn,
                              sG.number_of_nodes(), sG.number_of_edges())

                    id_pkl += 1
                    cursor.execute("""INSERT INTO subgraphs (
                                          id_pkl, 
                                          n_graphs, 
                                          graph6,
                                          k,
                                          k_partite,
                                          k_valences,
                                          nodes_valences,
                                          n_nodes, n_edges)
                                      VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""", (
                                          id_pkl,
                                          len(subgraphs[vn]),
                                          line_geng,
                                          len(p),
                                          str(p),
                                          str(vt),
                                          str(vn),
                                          sG.number_of_nodes(),
                                          sG.number_of_edges()))

                    with open(os.path.join(pkls_out, "{}.pkl".format(id_pkl)), "wb") as fn_pkls:
                        pickle.dump(root, fn_pkls)

            conn.commit()
    conn.close()
