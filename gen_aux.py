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
from databases import reformat_xml, create_connectivity_database


def get_from_hmdb(name, hmdb, out_dir):
    """
    Gets individual metabocards from HMDB.

    :param name: What to name the XML file.

    :param hmdb: The HMDB ID to be found.

    :param out_dir: The directory to save the XML file in.
    """

    urllib.request.urlretrieve("http://www.hmdb.ca/metabolites/" + hmdb + ".xml", os.path.join(out_dir, name + ".xml"))

    reformat_xml(os.path.join(out_dir, name + ".xml"))


def build_graph_isomorphism_database(db_out='../databases/k_graphs.sqlite', path=".."):
    #########################################
    # BUILD GRAPH ISOMORPHISM DATABASE
    #########################################
    # McKay, B.D. and Piperno, A., Practical Graph Isomorphism, II,
    # Journal of Symbolic Computation, 60 (2014), pp. 94-112, http://dx.doi.org/10.1016/j.jsc.2013.09.003
    # http://pallini.di.uniroma1.it/

    if sys.platform == "win32" or sys.platform == "win64":
        path_ri = os.path.join(path, "tools", "RI_win", "RI3.6-release", "ri36")

    elif sys.platform == "darwin":
        path_ri = os.path.join(path, "tools", "RI_mac", "RI3.6-release", "ri36")

    elif sys.platform == "linux2":
        if "bb" in "socket.gethostname":
            path_ri = os.path.join(path, "tools", "RI_unix", "RI3.6-release", "ri36")
        else:
            path_ri = os.path.join(path, "tools", "RI_bb", "RI3.6-release", "ri36")

    elif sys.platform == "linux":
        path_ri = os.path.join(path, "tools", "RI_unix", "RI3.6-release", "ri36")

    else:
        path_ri = os.path.join(path, "ri36")

    create_connectivity_database(db_out, path_ri=path_ri)
