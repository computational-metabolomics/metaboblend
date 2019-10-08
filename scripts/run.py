import os
import subprocess
import sys
import tempfile

sys.path.append(os.path.join("..", "metaboverse"))
from databases import update_substructure_database, create_isomorphism_database
from databases import SubstructureDb, ConnectivityDb
from build_structures import build


def main():

    #########################################
    # BUILD GRAPH ISOMORPHISM DATABASE
    #########################################
    # McKay, B.D. and Piperno, A., Practical Graph Isomorphism, II,
    # Journal of Symbolic Computation, 60 (2014), pp. 94-112, http://dx.doi.org/10.1016/j.jsc.2013.09.003
    # http://pallini.di.uniroma1.it/

    path = ".."
    if sys.platform == "win32" or sys.platform == "win64":
        path_geng = os.path.join(path, "tools", "nauty25r9_win", "geng")
        path_RI = os.path.join(path, "tools", "RI_win", "RI3.6-release", "ri36")

    elif sys.platform == "darwin":
        path_geng = os.path.join(path, "tools", "nauty25r9_mac", "geng")
        path_RI = os.path.join(path, "tools", "RI_mac", "RI3.6-release", "ri36")

    elif sys.platform == "linux2":
        if "bb" in "socket.gethostname":
            path_geng = os.path.join(path, "tools", "nauty25r9_unix", "geng")
            path_RI = os.path.join(path, "tools", "RI_unix", "RI3.6-release", "ri36")
        else:
            path_geng = os.path.join(path, "tools", "nauty25r9_bb", "geng")
            path_RI = os.path.join(path, "tools", "RI_bb", "RI3.6-release", "ri36")
    else:
        path_geng = os.path.join(path, "geng")
        path_RI = os.path.join(path, "ri36")

    path = os.path.join("..")
    db_out = os.path.join(path, "databases", "k_graphs.sqlite")
    pkls_out = os.path.join(path, "databases", "pkls")

    sizes = [1, 2]
    boxes = 3

    print("==================================")
    print(db_out)
    print(pkls_out)
    print(sizes, boxes)
    print("==================================")
    create_isomorphism_database(db_out, pkls_out, boxes, sizes, path_geng, path_RI)

    ##########################################
    # BUILD SUBSTRUCTURE DATABASE - HMDB
    ##########################################
    # 1-Methylhistidine - HMDB0000001
    # 2-Ketobutyric acid - HMDB0000005
    # 2-Hydroxybutyric - acid HMDB0000008
    # D-Glucose - HMDB0000122


    records = ["HMDB0000001", "HMDB0000005", "HMDB0000008", "HMDB0000122"]
    path_db = "../databases/substructures.sqlite"

    db = SubstructureDb(path_db, "", "")
    db.create_compound_database()
    db.close()

    for i, record in enumerate(records):
        HMDB_xml = os.path.join("input", record + ".xml")
        update_substructure_database(HMDB_xml, path_db, 2, 9)

    db = SubstructureDb(path_db, "", "")
    db.create_indexes()
    db.close()

    ##########################################
    # BUILD STRUCTURES
    ##########################################
    path_db_k_graphs = "../databases/k_graphs.sqlite"
    path_pkls = "../databases/pkls"

    path_out = "results/"
    db = SubstructureDb(path_db, path_pkls, path_db_k_graphs)

    # criteria substructures
    accuracy = 1
    heavy_atoms = range(2, 9)
    max_valence = 2

    # Select all HMDB compounds in database
    records = db.select_compounds()
    for record in records:

        print(record[0], str(record[11]))
        print(record[3:9], record[1])
        print("----------------------------------")

        fn_out = os.path.join(path_out, "{}.smi".format(record[0]))

        build(list(record[3:9]), record[1], db, fn_out, heavy_atoms, max_valence, accuracy)

        with open(fn_out) as smiles:
            temp_fn = tempfile.NamedTemporaryFile(mode="w", delete=False)
            for s in set(smiles.readlines()):
                temp_fn.write(s)

            # proc = subprocess.Popen(["obabel", "-ismi", os.path.join(path_out, "%s.smi" % (record[0])), "-osvg", "-O", os.path.join(path_out, "figures", "%s.svg" % (record[0]))], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            proc = subprocess.Popen(["obabel", "-ismi", os.path.join(path_out, temp_fn.name), "-osvg", "-O", os.path.join(path_out, "figures", "%s.svg" % (record[0]))], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #out, err = proc.communicate()
            #print out
            #print err
            temp_fn.seek(0)
            temp_fn.close()
    db.close()


if __name__ == "__main__":
    main()
