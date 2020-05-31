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
from databases import reformat_xml, update_substructure_database, filter_records, parse_xml, SubstructureDb, get_elements, calculate_exact_mass
from build_structures import build

sys.path.append(os.path.join("..", "functions"))
from gen_aux import get_uniq_subs

def hmdb_sdf_to_csv(path_in, path_out):
    """
    Converts the entire HMDB database, in sdf format, to a csv for processing by MetFrag.
    """

    with open(path_out, "w", encoding="utf-8", newline="") as csv_db_file:
        csv_db = csv.writer(csv_db_file)
        csv_db.writerow([
            "Identifier",
            "MonoisotopicMass",
            "MolecularFormula",
            "SMILES",
            "InChI",
            "InChIKey1",
            "InChIKey2",
            "InChIKey3",
            "Name",
            "InChIKey"
        ])

        for i, mol in enumerate(Chem.SDMolSupplier(path_in)):
            if mol is None:
                continue

            atom_check = [True for atom in mol.GetAtoms() if atom.GetSymbol() not in ["C", "H", "N", "O", "P", "S"]]
            if len(atom_check) > 0:
                continue

            Chem.SanitizeMol(mol)

            try:
                csv_db.writerow([
                    mol.GetProp("DATABASE_ID"),
                    mol.GetProp("EXACT_MASS"),
                    mol.GetProp("FORMULA"),
                    Chem.MolToSmiles(mol),
                    mol.GetProp("INCHI_IDENTIFIER"),
                    mol.GetProp("INCHI_KEY").split("-")[0],
                    mol.GetProp("INCHI_KEY").split("-")[1],
                    mol.GetProp("INCHI_KEY").split("-")[2],
                    mol.GetProp("GENERIC_NAME"),
                    mol.GetProp("INCHI_KEY")
                ])
            except KeyError:
                print(mol.GetProp("DATABASE_ID"))

    return path_out


def test_build(out_dir, mc, exact_mass, mol, hmdb_id, path_subs, path_k_graphs, path_pkls, heavy_atoms, max_valence,
               accuracy, fragment_masses, ppm, hydrogenation_allowance=2, max_atoms_available=2, max_n_substructures=3):
    """
    Altered build method - see msn_build.py for full description. Returns the recurrence dictionary for use in
    calculating the ROC graph.
    """

    mol_smi = Chem.MolToSmiles(mol, kekuleSmiles=True)
    pre_reccurence = {}

    for fragment_mass in fragment_masses:
        smi_out = os.path.join(out_dir, "{}_".format(hmdb_id) + str(round(fragment_mass, 4)) + ".smi")
        open(smi_out, "w").close()

        for j in range(0 - hydrogenation_allowance, hydrogenation_allowance + 1):
            hydrogenated_fragment_mass = fragment_mass + (j * 1.007825)
            build(mc, exact_mass, smi_out, heavy_atoms, max_valence, accuracy, max_atoms_available, max_n_substructures,
                  path_k_graphs, path_pkls, path_subs, hydrogenated_fragment_mass, ppm, out_mode="a")

        get_uniq_subs(smi_out, ignore_substructures=True)
        with open(smi_out, mode="r") as smis:
            for line in smis:
                if len(line) > 0:
                    try:
                        pre_reccurence[line.strip()] += 1
                    except KeyError:
                        pre_reccurence[line.strip()] = 1

    with open(os.path.join(out_dir, "structure_ranks.csv"), "w", newline="") as ranks_out:
        ranks_csv = csv.writer(ranks_out)
        ranks_csv.writerow(["kekule_smiles", "occurence"])

        num_recurrence = {}
        num_struct = 0
        for smi in pre_reccurence.keys():
            ranks_csv.writerow([smi, pre_reccurence[smi]])
            num_struct += 1
            try:
                num_recurrence[pre_reccurence[smi]] += 1

            except KeyError:
                num_recurrence[pre_reccurence[smi]] = 1
    try:
        recurrence = str(pre_reccurence[mol_smi])

    except KeyError:
        recurrence = "0"
        better_candidates = 0
        for num in num_recurrence.keys():
            better_candidates += num_recurrence[num]

    else:
        better_candidates = 0
        for num in num_recurrence.keys():
            if num >= pre_reccurence[mol_smi]:
                better_candidates += num_recurrence[num]

    if len(pre_reccurence.values()) > 0:
        max_recurrence = str(max(pre_reccurence.values()))
    else:
        max_recurrence = str(0)

    return pre_reccurence