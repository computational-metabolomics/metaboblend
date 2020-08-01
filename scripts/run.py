from metaboblend.databases import *
from metaboblend.build_structures import *


def build_graph_isomorphism_database(sizes=[1, 2], boxes=3, db_out='../databases/k_graphs.sqlite'):
    #########################################
    # BUILD GRAPH ISOMORPHISM DATABASE
    #########################################
    # McKay, B.D. and Piperno, A., Practical Graph Isomorphism, II,
    # Journal of Symbolic Computation, 60 (2014), pp. 94-112, http://dx.doi.org/10.1016/j.jsc.2013.09.003
    # http://pallini.di.uniroma1.it/

    if sys.platform == "win32" or sys.platform == "win64":
        path_ri = os.path.join("..", "tools", "RI_win", "RI3.6-release", "ri36")

    elif sys.platform == "darwin":
        path_ri = os.path.join("..", "tools", "RI_mac", "RI3.6-release", "ri36")

    elif sys.platform == "linux2":
        if "bb" in "socket.gethostname":
            path_ri = os.path.join("..", "tools", "RI_unix", "RI3.6-release", "ri36")
        else:
            path_ri = os.path.join("..", "tools", "RI_bb", "RI3.6-release", "ri36")

    elif sys.platform == "linux":
        path_ri = os.path.join("..", "tools", "RI_unix", "RI3.6-release", "ri36")

    else:
        path_ri = os.path.join("..", "ri36")

    print("==================================")
    print(db_out)
    print(sizes, boxes)
    print("==================================")

    create_isomorphism_database(db_out, boxes, sizes, path_ri)


def build_structures(max_n_substructures=3, path_db_k_graphs="../databases/k_graphs.sqlite",
                     path_pkls="../databases/pkls", path_db="../databases/substructures.sqlite", path_out="results/"):

    db = SubstructureDb(path_db, path_pkls, path_db_k_graphs)

    # Select all HMDB compounds in database
    records = db.select_compounds()
    for record in records:
        print(record)
        print("----------------------------------")

        fn_out = os.path.join(path_out, "{}.smi".format(record[0]))

        # build & write structures
        build(list(record[3:9]), record[1], max_n_substructures, fn_out, path_connectivity_db=path_db_k_graphs,
              path_substructure_db=path_db, out_mode="w", ppm=None, processes=None, table_name="substructures",
              prescribed_mass=None)

    db.close()


if __name__ == "__main__":
    build_graph_isomorphism_database()

    hmdbs = [os.path.join("input", h + ".xml") for h in ["HMDB0000001", "HMDB0000005", "HMDB0000008", "HMDB0000122"]]
    create_substructure_database(hmdbs, "../databases/substructures.sqlite", 2, 9)

    build_structures()
