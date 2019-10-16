import io
import os
from io import BytesIO
import subprocess
import pickle
import pprint
import sqlite3
import sys
import tempfile
from collections import OrderedDict
import xml.etree.ElementTree as etree
import networkx as nx
from rdkit import Chem
from auxiliary import calculate_complete_multipartite_graphs, graph_to_ri, graph_info

sqlite3.register_converter("PICKLE", pickle.loads)


def parse_xml(source, encoding="utf8"):

    with io.open(source, "r", encoding=encoding) as inp:
        record_out = OrderedDict()

        xmldec = inp.readline()
        xmldec2 = inp.readline()

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

    def __init__(self, db):
        self.db = db


class SubstructureDb:

    def __init__(self, db, path_pkls, db2=None):
        self.db = db
        self.db2 = db2
        self.path_pkls = path_pkls

        self.conn = sqlite3.connect(self.db)
        self.cursor = self.conn.cursor()

        if self.db2 is not None:
            self.cursor.execute("""ATTACH DATABASE '%s' as 'graphs';""" % self.db2)

    def select_compounds(self, cpds=[]):
        if len(cpds) > 0:
            sql = " WHERE HMDBID in ('%s')" % (", ".join(map(str, cpds)))
        else:
            sql = ""

        self.cursor.execute("""select distinct HMDBID, exact_mass, formula, C, H, N, O, P, S, SMILES,
                            SMILES_RDKIT, SMILES_RDKIT_KEK from compounds%s""" % sql)
        return self.cursor.fetchall()

    def select_mass_values(self, accuracy, heavy_atoms, max_valence, masses, db):
        mass_values = []
        filter_mass = ""
        if type(masses) == list:
            if len(masses) > 0:
                filter_mass = " AND exact_mass__1 in ({})".format(",".join(map(str, masses)))

        self.cursor.execute("""SELECT DISTINCT exact_mass__{}
                                   FROM substructures 
                               WHERE valence <= {}
                                   AND heavy_atoms IN ({}){}
                            """.format(accuracy, max_valence, ",".join(map(str, heavy_atoms)), filter_mass))

        records = self.cursor.fetchall()
        for record in records:
            mass_values.append(record[0])
        mass_values.sort()
        return mass_values

    def select_ecs(self, exact_mass, heavy_atoms, accuracy):
        self.cursor.execute("""SELECT DISTINCT 
                                   C, 
                                   H, 
                                   N, 
                                   O, 
                                   P, 
                                   S 
                               FROM substructures 
                               WHERE heavy_atoms in ({})
                               AND exact_mass__{} = {}
                            """.format(",".join(map(str, heavy_atoms)), accuracy, exact_mass))
        return self.cursor.fetchall()

    def paths(self, tree, cur=()):
        if tree == {}:
            yield cur
        else:
            for n, s in tree.items():
                for path in self.paths(s, cur+(n,)):
                    yield path

    def isomorphism_graphs(self, id_pkl):
        with open(os.path.join(self.path_pkls, "{}.pkl".format(id_pkl)), 'rb') as pickle_file:
            nGcomplete = pickle.load(pickle_file)
        for p in self.paths(nGcomplete):
            yield p

    def k_configs(self):
        self.cursor.execute("""SELECT id_pkl, nodes_valences 
                               FROM subgraphs""")
        records = self.cursor.fetchall()
        configs = {}
        for record in records:
            configs[str(record[1])] = record[0]
        return configs

    def select_sub_structures(self, l_atoms):

        subsets = []
        for i in range(len(l_atoms)):

            self.cursor.execute("""SELECT DISTINCT lib 
                                   FROM substructures
                                   WHERE C = {} 
                                   AND H = {} 
                                   AND N = {} 
                                   AND O = {}
                                   AND P = {}
                                   AND S = {}
                                """.format(l_atoms[i][0], l_atoms[i][1], l_atoms[i][2], l_atoms[i][3], l_atoms[i][4], l_atoms[i][5]))
            records = self.cursor.fetchall()
            if len(records) == 0:
                return []
            ss = [pickle.loads(record[0]) for record in records]
            subsets.append(ss)

        return subsets

    def create_compound_database(self):
        self.cursor.execute('DROP TABLE IF EXISTS compounds')
        self.cursor.execute('DROP TABLE IF EXISTS substructures')
        self.cursor.execute('DROP TABLE IF EXISTS hmdbid_substructures')

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
                              smiles_rdkit_kek,
                              PRIMARY KEY (hmdbid, smiles_rdkit_kek))""")

    def create_indexes(self):

        self.cursor.execute("""DROP INDEX IF EXISTS heavy_atoms__Valence__mass__1__idx""")
        self.cursor.execute("""DROP INDEX IF EXISTS heavy_atoms__Valence__mass__0_1__idx""")
        self.cursor.execute("""DROP INDEX IF EXISTS heavy_atoms__Valence__mass__0_01__idx""")
        self.cursor.execute("""DROP INDEX IF EXISTS heavy_atoms__Valence__mass__0_001__idx""")
        self.cursor.execute("""DROP INDEX IF EXISTS heavy_atoms__Valence__mass__0_0001__idx""")
        self.cursor.execute("""DROP INDEX IF EXISTS atoms__Valence__idx""")

        self.cursor.execute("""CREATE INDEX heavy_atoms__Valence__mass__1__idx 
                               ON substructures (heavy_atoms, valence, valence_atoms, exact_mass__1);""")
        self.cursor.execute("""CREATE INDEX heavy_atoms__Valence__mass__0_1__idx 
                               ON substructures (heavy_atoms, valence, valence_atoms, exact_mass__0_1);""")
        self.cursor.execute("""CREATE INDEX heavy_atoms__Valence__mass__0_01__idx 
                               ON substructures (heavy_atoms, valence, valence_atoms, exact_mass__0_01);""")
        self.cursor.execute("""CREATE INDEX heavy_atoms__Valence__mass__0_001__idx 
                               ON substructures (heavy_atoms, valence, valence_atoms, exact_mass__0_001);""")
        self.cursor.execute("""CREATE INDEX heavy_atoms__Valence__mass__0_0001__idx
                               ON substructures (heavy_atoms, valence, valence_atoms, exact_mass__0_0001);""")
        self.cursor.execute("""CREATE INDEX atoms__Valence__idx 
                               ON substructures (C, H, N, O, P, S, valence, valence_atoms);""")

    def close(self):
        self.conn.close()


def get_substructure(mol, idxs_edges_subgraph, debug=False):

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

    # Returns the type of the bond as a double (i.e. 1.0 for SINGLE, 1.5 for AROMATIC, 2.0 for DOUBLE)

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

                if atom_n.GetIdx() not in degree_atoms:
                    degree_atoms[atom_n.GetIdx()] = 1
                else:
                    degree_atoms[atom_n.GetIdx()] += 1

    bond_types = {}

    for b in mol_out.GetBonds():
        if debug:
            print(b.GetBondTypeAsDouble())
            print(b.GetBondType())

        print(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), mol_out.GetAtomWithIdx(b.GetBeginAtomIdx()).GetSymbol(), mol_out.GetAtomWithIdx(b.GetEndAtomIdx()).GetSymbol())
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
        Chem.rdmolops.Kekulize(mol_out)
    except:
        return None

    return {"smiles": Chem.MolToSmiles(mol_out, kekuleSmiles=True), # REORDERED ATOM INDEXES
            "mol": mol_out,
            "bond_types": bond_types,
            "degree_atoms": degree_atoms,
            "valence": sum(degree_atoms.values()),
            "atoms_available": len(degree_atoms.keys()),
            "dummies": dummies}


def get_elements(mol, elements=None):
    if not elements:
        elements = {"C": 0, "H": 0, "N": 0, "O": 0, "P": 0, "S": 0, "*": 0}
    mol = Chem.AddHs(mol)
    for atom in mol.GetAtoms():
        elements[atom.GetSymbol()] += 1
    return elements


def calculate_exact_mass(mol, exact_mass_elements=None):
    if not exact_mass_elements:
        exact_mass_elements = {"C": 12.0, "H": 1.007825, "N": 14.003074, "O": 15.994915, "P": 30.973763, "S": 31.972072, "*": -1.007825}
    exact_mass = 0.0
    mol = Chem.AddHs(mol)
    for atom in mol.GetAtoms():
        atomSymbol = atom.GetSymbol()
        if atomSymbol != "*":
            exact_mass += exact_mass_elements[atomSymbol]
    return exact_mass


def filter_records(records, hmdb_cpd=None):

    for record in records:

        if "smiles" in record:

            mol = Chem.MolFromSmiles(record['smiles'])

            if mol is None:
                continue

            if mol.GetNumHeavyAtoms() < 4:
                continue

            atom_check = [True for atom in mol.GetAtoms() if atom.GetSymbol() not in ["C", "H", "N", "O", "P", "S"]]
            if len(atom_check) > 0:
                continue

            smiles = Chem.rdmolfiles.MolToSmiles(mol, kekuleSmiles=True)
            Chem.rdmolops.Kekulize(mol)
            smiles_rdkit_kek = Chem.rdmolfiles.MolToSmiles(mol, kekuleSmiles=True)

            if "+" in smiles_rdkit_kek or "-" in smiles_rdkit_kek or "+" in smiles or "-" in smiles:
                # print record['HMDB_ID'], record['smiles'], "+/-"
                continue

            print("%s\t%s" % (record['accession'], record['monisotopic_molecular_weight']))

            els = get_elements(mol)
            exact_mass = calculate_exact_mass(mol)

            record_dict = {'HMDB_ID': record['accession'],
                           'formula': record["chemical_formula"],
                           'exact_mass': round(exact_mass, 6),
                           'smiles': record['smiles'],
                           'smiles_rdkit': smiles,
                           'smiles_rdkit_kek': smiles_rdkit_kek,
                           'C': els['C'],
                           'H': els['H'],
                           'N': els['N'],
                           'O': els['O'],
                           'P': els['P'],
                           'S': els['S'],
                           'mol': mol}

            yield record_dict


def update_substructure_database(fn_hmdb, fn_db, n_min, n_max, hmdb_cpd=None):

    conn = sqlite3.connect(fn_db)
    cursor = conn.cursor()

    records = parse_xml(fn_hmdb)

    for record_dict in filter_records(records, hmdb_cpd):

        cursor.execute("""INSERT OR IGNORE INTO compounds (
                              hmdbid, 
                              exact_mass, 
                              formula, 
                              C, H, N, O, P, S, 
                              smiles, 
                              smiles_rdkit, 
                              smiles_rdkit_kek)
                          values (
                              :HMDB_ID, 
                              :exact_mass,
                              :formula, 
                              :C, :H, :N, :O, :P, :S, 
                              :smiles, 
                              :smiles_rdkit, 
                              :smiles_rdkit_kek)""", record_dict)

        # Returns a tuple of 2-tuples with bond IDs

        for sgs in Chem.rdmolops.FindAllSubgraphsOfLengthMToN(record_dict["mol"], n_min, n_max):
            for edge_idxs in sgs:
                lib = get_substructure(record_dict["mol"], edge_idxs)
                if lib is None:
                    continue

                smiles_rdkit_kek = Chem.rdmolfiles.MolToSmiles(lib["mol"], kekuleSmiles=True)

                exact_mass = calculate_exact_mass(lib["mol"])
                els = get_elements(lib["mol"])

                pkl_lib = pickle.dumps(lib)
                sub_smi_dict = {'smiles': smiles_rdkit_kek,
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
                                      :atoms_available,:lib)
                               """, sub_smi_dict)

                cursor.execute("""INSERT OR IGNORE INTO hmdbid_substructures (
                                      hmdbid, 
                                      smiles_rdkit_kek) 
                                  VALUES ("%s", "%s")
                               """ % (record_dict['HMDB_ID'], smiles_rdkit_kek))
    conn.commit()
    conn.close()


def create_isomorphism_database(db_out, pkls_out, boxes, sizes, path_geng=None, path_RI=None):

    conn = sqlite3.connect(db_out)
    cursor = conn.cursor()

    cursor.execute('''DROP TABLE IF EXISTS subgraphs''')
    cursor.execute('''CREATE TABLE subgraphs (
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
                   );''')
    conn.commit()

    id_pkl = 0

    for G, p in calculate_complete_multipartite_graphs(sizes, boxes):

        print([path_geng, str(G.number_of_nodes()), "-d1", "-D2", "-q"])
        proc = subprocess.Popen([path_geng, str(len(G.nodes)), "-d1", "-D2", "-q"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        geng_out, err = proc.communicate()

        proc.stdout.close()
        proc.stderr.close()

        for i, line_geng in enumerate(geng_out.split()):

            print(line_geng)

            sG = nx.read_graph6(BytesIO(line_geng))

            k_gfu = tempfile.NamedTemporaryFile(mode="w", delete=False)
            k_gfu.write(graph_to_ri(G, "k_graph"))
            k_gfu.seek(0)

            s_gfu = tempfile.NamedTemporaryFile(mode="w", delete=False)
            s_gfu.write(graph_to_ri(sG, "subgraph"))
            s_gfu.seek(0)

            proc = subprocess.Popen([path_RI, "mono", "geu", k_gfu.name, s_gfu.name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            RI_out, err = proc.communicate()

            k_gfu.close()
            s_gfu.close()

            mappings = []
            subgraphs = {}

            for line in RI_out.decode("utf-8").splitlines():
                if line[0] == "{":
                    mappings.append(eval(line))

                if len(mappings) == 20000:
                    gi = graph_info(p, sG, mappings,)

                    for vn in gi[0]:

                        if vn not in subgraphs:
                            subgraphs[vn] = gi[0][vn]
                            # print vn, result[0][vn], result[1][0], result[1][1], len(result[1][1])
                        else:

                            before = len(subgraphs[vn])
                            for es in gi[0][vn]:
                                if es not in subgraphs[vn]:
                                    subgraphs[vn].append(es)
                                    # print vn, es, result[1][0], result[1][1], len(result[1][1])
                            after = len(subgraphs[vn])
                            print(before, after)

                    mappings = []

            if len(mappings) > 0:
                gi = graph_info(p, sG, mappings, )
                #job = job_server.submit(graphInfo, (p, sG, mappings, ), (valences,), modules=(), globals=globals())
                #jobs.append(job)

                for vn in gi[0]:

                    if vn not in subgraphs:
                        subgraphs[vn] = gi[0][vn]
                        #print vn, result[0][vn], result[1][0], result[1][1], len(result[1][1])
                    else:

                        before = len(subgraphs[vn])
                        for es in gi[0][vn]:
                            if es not in subgraphs[vn]:
                                subgraphs[vn].append(es)
                                #print vn, es, result[1][0], result[1][1], len(result[1][1])
                        after = len(subgraphs[vn])
                        print(before, after)

            if len(subgraphs) > 0:

                for vn in subgraphs:

                    root = {}
                    for fr in subgraphs[vn]:
                        parent = root
                        for e in fr:
                            parent = parent.setdefault(e, {})

                    vt = tuple([sum(v) for v in eval(vn)])
                    print("INSERT:", i, line_geng.decode("utf-8"), len(subgraphs[vn]), len(p), str(p), vt, vn, sG.number_of_nodes(), sG.number_of_edges())

                    id_pkl += 1
                    cursor.execute('''INSERT INTO subgraphs (id_pkl, 
                                      n_graphs, 
                                      graph6,
                                      k,
                                      k_partite,
                                      k_valences,
                                      nodes_valences,
                                      n_nodes, n_edges) 
                                      values (?, ?, ?, ?, ?, ?, ?, ?, ?)''', (
                                      id_pkl,
                                      len(subgraphs[vn]),
                                      line_geng,
                                      len(p),
                                      str(p),
                                      str(vt),
                                      str(vn),
                                      sG.number_of_nodes(),
                                      sG.number_of_edges()))
                    pickle.dump(root, open(os.path.join(pkls_out, "{}.pkl".format(id_pkl)), "wb"))
            conn.commit()
    conn.close()
