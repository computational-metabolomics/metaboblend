import os
import sys
import urllib.request
import csv
import tempfile
import networkx as nx
from rdkit import Chem
from shutil import rmtree
import pickle

sys.path.append(os.path.join("..", "..", "..", "metaboblend", "metaboblend"))
from databases import reformat_xml, update_substructure_database, filter_records, parse_xml, SubstructureDb, get_elements, calculate_exact_mass, create_isomorphism_database
from build_structures import build


def get_from_hmdb(name, hmdb, out_dir):
    """
    Gets individual metabocards from HMDB.

    :param name: What to name the XML file.

    :param hmdb: The HMDB ID to be found.

    :param out_dir: The directory to save the XML file in.
    """

    urllib.request.urlretrieve("http://www.hmdb.ca/metabolites/" + hmdb + ".xml", os.path.join(out_dir, name + ".xml"))

    reformat_xml(os.path.join(out_dir, name + ".xml"))


def build_substructure_database(records, path_input, path_db="../databases/substructures.sqlite", n_min=2,
                                n_max=9, method="exhaustive", max_atoms_available=None, max_valence=None):
    """
    Build a substructure database from a series of compounds.

    :param records: HMDB IDs to be used to generate the substructure database.

    :param path_input: The path of HMDB XML files.

    :param path_db: The path at which to generate the SQLite database.

    :param n_min: The minimum number of bonds to be included in substructures.

    :param n_max: The maximum number of bonds to be included in substructures.

    :param method: The fragmentation method to be used.
    """

    db = SubstructureDb(path_db, "", "")
    db.create_compound_database()
    db.close()

    for i, record in enumerate(records):
        update_substructure_database(os.path.join(path_input, record + ".xml"), path_db, n_min, n_max, method=method,
                                     max_atoms_available=max_atoms_available, max_valence=max_valence)

    db = SubstructureDb(path_db, "", "")
    db.create_indexes()
    db.close()


def get_uniq_subs(smi_out, ignore_substructures=False):
    """
    Filter out multiply generated substructures.

    :param smi_out: The file of smiles to be filtered.

    :param ignore_substructures: If True, get unique structures, if False, get unique combinations of structures and
        the substructures they were generated from.
    """
    temp_fn = tempfile.NamedTemporaryFile(mode="w")
    seen_lines = set()

    with open(smi_out, "r") as structures:
        for line in structures:
            if ignore_substructures:
                if line.split()[0] not in seen_lines:
                    seen_lines.add(line.split()[0] + "\n")
            else:
                if line not in seen_lines:
                    seen_lines.add(line)

    temp_fn.close()

    with open(smi_out, "w") as unique_structures:
        for line in seen_lines:
            unique_structures.writelines(line)


def subset_substructures(hmdb_ids, in_db, out_db, substructures_only=True, subset=True):
    """
    Take the substructures of a single, or series of compounds, from one database, and create a new substructure
    database that only includes them.

    :param hmdb_ids: List of HMDB IDs.

    :param in_db: The database to extract substructures from.

    :param out_db: The new database to be generated.

    :param substructures_only: Only extract the substructure database if True, else extract all major tables.

    :param subset: If True, subset the substructures by the HMDBIDs, else simply makes a copy of the database.
    """

    original_substructures = SubstructureDb(in_db, "")

    original_substructures.cursor.execute("attach database '%s' as subset" % out_db)

    original_substructures.cursor.execute("drop table if exists subset.substructures")
    original_substructures.cursor.execute("drop table if exists subset.compounds")
    original_substructures.cursor.execute("drop table if exists subset.hmdbid_substructures")

    if subset:
        original_substructures.cursor.execute("""create table subset.substructures as 
                                                 select * from substructures 
                                                 where smiles in (
                                                    select smiles from hmdbid_substructures
                                                    where hmdbid in ({seq})
                                                 )""".format(seq=','.join(['?']*len(hmdb_ids))), hmdb_ids)
    else:
        original_substructures.cursor.execute("""create table subset.substructures as 
                                                 select * from substructures""")

    original_substructures.cursor.execute("""select * from substructures 
                                                 where smiles in (
                                                    select smiles from hmdbid_substructures
                                                    where hmdbid in ({seq})
                                                 )""".format(seq=','.join(['?'] * len(hmdb_ids))), hmdb_ids)

    if not substructures_only:
        original_substructures.cursor.execute("""create table subset.compounds as 
                                                 select * from compounds
                                                 where hmdbid in ({seq}
                                                 )""".format(seq=','.join(['?']*len(hmdb_ids))), hmdb_ids)

        original_substructures.cursor.execute("""create table subset.hmdbid_substructures as 
                                                 select * from hmdbid_substructures
                                                 where hmdbid in ({seq}
                                                 )""".format(seq=','.join(['?']*len(hmdb_ids))), hmdb_ids)


def build_graph_isomorphism_database(sizes=[1, 2], boxes=3,
                                     db_out='../databases/k_graphs.sqlite',
                                     pkls_out='../databases/pkls', path="..", debug=False):
    #########################################
    # BUILD GRAPH ISOMORPHISM DATABASE
    #########################################
    # McKay, B.D. and Piperno, A., Practical Graph Isomorphism, II,
    # Journal of Symbolic Computation, 60 (2014), pp. 94-112, http://dx.doi.org/10.1016/j.jsc.2013.09.003
    # http://pallini.di.uniroma1.it/

    if sys.platform == "win32" or sys.platform == "win64":
        path_geng = os.path.join(path, "tools", "nauty25r9_win", "geng")
        path_ri = os.path.join(path, "tools", "RI_win", "RI3.6-release", "ri36")

    elif sys.platform == "darwin":
        path_geng = os.path.join(path, "tools", "nauty25r9_mac", "geng")
        path_ri = os.path.join(path, "tools", "RI_mac", "RI3.6-release", "ri36")

    elif sys.platform == "linux2":
        if "bb" in "socket.gethostname":
            path_geng = os.path.join(path, "tools", "nauty25r9_unix", "geng")
            path_ri = os.path.join(path, "tools", "RI_unix", "RI3.6-release", "ri36")
        else:
            path_geng = os.path.join(path, "tools", "nauty25r9_bb", "geng")
            path_ri = os.path.join(path, "tools", "RI_bb", "RI3.6-release", "ri36")

    elif sys.platform == "linux":
        path_geng = os.path.join(path, "tools", "nauty25r9_unix", "geng")
        path_ri = os.path.join(path, "tools", "RI_unix", "RI3.6-release", "ri36")

    else:
        path_geng = os.path.join(path, "geng")
        path_ri = os.path.join(path, "ri36")

    print("==================================")
    print(db_out)
    print(pkls_out)
    print(sizes, boxes)
    print("==================================")
    create_isomorphism_database(db_out, pkls_out, boxes, sizes, path_geng, path_ri, debug=debug)
