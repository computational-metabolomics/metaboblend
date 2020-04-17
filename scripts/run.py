import os
import subprocess
import sys
import tempfile
from metaboverse.databases import *


def build_graph_isomorphism_database(sizes=[1, 2], boxes=3,
                                     db_out='../databases/k_graphs.sqlite',
                                     pkls_out='../databases/pkls', path=".."):
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
    create_isomorphism_database(db_out, pkls_out, boxes, sizes, path_geng, path_ri)


def build_substructure_database(records, path_input, path_db="../databases/substructures.sqlite", n_min=2,
                                n_max=9, method="exhaustive"):
    db = SubstructureDb(path_db, "", "")
    db.create_compound_database()
    db.close()

    for i, record in enumerate(records):
        update_substructure_database(os.path.join(path_input, record + ".xml"), path_db, n_min, n_max, method=method)

    db = SubstructureDb(path_db, "", "")
    db.create_indexes()
    db.close()


def build_structures(accuracy=1, heavy_atoms=range(2, 9), max_valence=4, max_atoms_available=2, max_n_substructures=3,
                     path_db_k_graphs="../databases/k_graphs.sqlite",
                     path_pkls="../databases/pkls",
                     path_db="../databases/substructures.sqlite",
                     path_out="results/", debug=False):
    db = SubstructureDb(path_db, path_pkls, path_db_k_graphs)

    # Select all HMDB compounds in database
    records = db.select_compounds()
    for record in records:
        print(record[0], str(record[11]))
        print(record[3:9], record[1])
        print("----------------------------------")

        fn_out = os.path.join(path_out, "{}.smi".format(record[0]))

        # build & write structures
        build(list(record[3:9]), record[1], db, fn_out, heavy_atoms, max_valence, accuracy, max_atoms_available,
              max_n_substructures, debug=debug)

        # write figures to svg files
        with open(fn_out) as smiles:
            temp_fn = tempfile.NamedTemporaryFile(mode="w", delete=False)
            for s in set(smiles.readlines()):
                temp_fn.write(s)

            subprocess.Popen(["obabel", "-ismi", os.path.join(path_out, temp_fn.name), "-osvg", "-O",
                              os.path.join(path_out, "figures", "%s.svg" % (record[0]))],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            temp_fn.seek(0)
            temp_fn.close()
    db.close()


if __name__ == "__main__":
    build_graph_isomorphism_database()
    build_substructure_database(["HMDB0000001", "HMDB0000005", "HMDB0000008", "HMDB0000122"], "input")
    build_structures()
