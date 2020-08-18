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
