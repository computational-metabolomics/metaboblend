import os
import sys
import urllib.request
import csv
import tempfile
import networkx as nx
from statistics import median
from rdkit import Chem
from shutil import rmtree
import pickle
from msp2db.db import create_db
from msp2db.parse import LibraryData
import sqlite3

metaboverse_path = os.path.join("..", "..", "..", "metaboverse")

sys.path.append(os.path.join(metaboverse_path, "metaboverse"))
from databases import update_substructure_database, filter_records, parse_xml, SubstructureDb, get_elements, calculate_exact_mass
from build_structures import build


def parse_testing_data(csv_path, hmdb_path):
    with open(csv_path) as csv_file:
        csv_data = csv.reader(csv_file)

        data_categories = {}
        for line in csv_data:
            if line[0] + ".xml" not in os.listdir(hmdb_path):
                get_from_hmdb(line[0], line[0], hmdb_path)

            # setup ms data category
            if line[4] not in data_categories:
                data_categories[line[4]] = {}

            # setup compound dict
            if line[0] not in data_categories[line[4]]:
                data_categories[line[4]][line[0]] = {"name": line[1], "precursor_ion_mass": float(line[3]), "peaks": []}

                for record_dict in filter_records(parse_xml(os.path.join(hmdb_path, line[0] + ".xml"))):
                    assert record_dict["HMDB_ID"] == line[0]

                    data_categories[line[4]][line[0]]["mc"] = [record_dict["C"], record_dict["H"], record_dict["N"],
                                                               record_dict["O"], record_dict["P"], record_dict["S"]]
                    data_categories[line[4]][line[0]]["exact_mass"] = record_dict["exact_mass"]
                    mol = Chem.MolFromSmiles(Chem.MolToSmiles(record_dict["mol"]))
                    Chem.SanitizeMol(mol)
                    data_categories[line[4]][line[0]]["mol"] = mol
                    data_categories[line[4]][line[0]]["smiles"] = Chem.MolToSmiles(mol)

            assert data_categories[line[4]][line[0]]["exact_mass"] is not None
            data_categories[line[4]][line[0]]["peaks"].append(float(line[2]))

    return data_categories


def add_small_substructures(path_subs):
    substructures = SubstructureDb(path_subs, "")
    small_mols = ["C*", "*C*", "O*", "*O*", "N*", "*N*"]
    for smi in small_mols:
        mol = Chem.MolFromSmiles(smi)
        Chem.SanitizeMol(mol)

        degree_atoms = {}
        dummies = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "*"]

        for atom in mol.GetAtoms():

            if atom.GetIdx() in dummies:

                for atom_n in atom.GetNeighbors():

                    if atom_n.GetSymbol() == "*":
                        continue  # do not count dummies for valence calculations
                    elif atom_n.GetIdx() not in degree_atoms:
                        degree_atoms[atom_n.GetIdx()] = 1
                    else:
                        degree_atoms[atom_n.GetIdx()] += 1

        bond_types = {}

        for b in mol.GetBonds():
            if mol.GetAtomWithIdx(b.GetBeginAtomIdx()).GetSymbol() == "*":
                if b.GetEndAtomIdx() not in bond_types:
                    bond_types[b.GetEndAtomIdx()] = [b.GetBondTypeAsDouble()]
                else:
                    bond_types[b.GetEndAtomIdx()].append(b.GetBondTypeAsDouble())

            elif mol.GetAtomWithIdx(b.GetEndAtomIdx()).GetSymbol() == "*":
                if b.GetBeginAtomIdx() not in bond_types:
                    bond_types[b.GetBeginAtomIdx()] = [b.GetBondTypeAsDouble()]
                else:
                    bond_types[b.GetBeginAtomIdx()].append(b.GetBondTypeAsDouble())

        lib = {"smiles": Chem.MolToSmiles(mol),
            "mol": mol,
            "bond_types": bond_types,
            "degree_atoms": degree_atoms,
            "valence": sum(degree_atoms.values()),
            "atoms_available": len(degree_atoms.keys()),
            "dummies": dummies}

        smiles_rdkit_kek = Chem.rdmolfiles.MolToSmiles(lib["mol"])

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

        substructures.cursor.execute("""INSERT OR IGNORE INTO substructures (
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

    substructures.conn.commit()
    substructures.close()


class MspDatabase:
    def __init__(self, path_db, path_msp=None, schema="mona"):
        self.path_db = path_db
        self.schema = schema

        if path_msp is not None:
            print("Generating DB")
            self.lib_data = self.generate_db(path_msp)
        else:
            assert os.path.exists(path_db)
            self.lib_data = None
            print("Using existing DB")

        self.conn = sqlite3.connect(self.path_db)
        self.cursor = self.conn.cursor()

    def generate_db(self, path_msp):
        create_db(file_pth=self.path_db)

        return LibraryData(msp_pth=path_msp, db_pth=self.path_db, db_type='sqlite', schema=self.schema)

    def get_fragments(self, precursor_type="[M+H]+", ms_level=2, max_mass=200, min_mass=100, snr=2.0):

        self.cursor.execute("""select distinct inchikey_id, name, exact_mass, smiles
                               from metab_compound where exact_mass > {0} and exact_mass < {1}
                            """.format(str(min_mass), str(max_mass)))

        compounds = self.cursor.fetchall()

        for compound in compounds:
            self.cursor.execute("""select id, accession, precursor_mz, polarity, precursor_type 
                                   from library_spectra_meta 
                                   where inchikey_id = '{0}' and ms_level = {1} and precursor_type = '{2}'
                                   """.format(compound[0], str(ms_level), precursor_type))

            spectra_meta = self.cursor.fetchone()
            if spectra_meta is None:
                continue
            elif len(spectra_meta) == 0:
                continue

            mz = []
            intensity = []

            self.cursor.execute("""select mz, i from library_spectra 
                                   where library_spectra_meta_id = %s""" % spectra_meta[0])

            for spectra in self.cursor.fetchall():
                mz.append(spectra[0])
                intensity.append(spectra[1])

            # inchikey_id, name, exact_mass, smiles, accession, precursor_mz, mzs, itensities
            # med_snr = median(intensity) * snr
            # indices = [i for i in range(len(mz)) if intensity[i] > med_snr]
            indices = sorted(range(len(intensity)), key=lambda i: intensity[i])[-15:]  # top 15 results
            mz, intensity = [m for i, m in enumerate(mz) if i in indices], \
                            [inten for i, inten in enumerate(intensity) if i in indices]

            if len(intensity) == 0:
                continue

            yield list(compound) + list(spectra_meta[1:3]) + [mz] + [intensity]

    def close(self):
        self.conn.close()


def parse_msp_testing_data(paths_msp_db, names_msp, path_hmdb_ids):
    with open(path_hmdb_ids) as hmdb_ids:
        hmdb_ids_csv = csv.reader(hmdb_ids)
        
        hmdb_dict = {}
        for row in hmdb_ids_csv:
            hmdb_dict[row[1]] = row[4]
            
    seen_hmdbs = set()
    data_categories = {}
    for path_msp_db, name_msp in zip(paths_msp_db, names_msp):
        data_categories[name_msp] = {}

        msp_db = MspDatabase(path_msp_db)

        # 0            1     2           3       4          5             6    7
        # inchikey_id, name, exact_mass, smiles, accession, precursor_mz, mzs, itensities
        msp_data = msp_db.get_fragments(precursor_type="[M+H]+")

        for spectra in msp_data:
            if spectra[0] in data_categories[name_msp]:
                continue

            mol = Chem.MolFromSmiles(spectra[3])
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

            if "+" in Chem.MolToSmiles(mol) or "-" in Chem.MolToSmiles(mol):
                continue

            try:
                hmdb_id = hmdb_dict[spectra[0]]
            except KeyError:
                continue

            if hmdb_id in seen_hmdbs:
                continue
            else:
                seen_hmdbs.add(hmdb_id)

            data_categories[name_msp][hmdb_id] = {"name": spectra[1], "inchikey_id": spectra[0],
                                                      "precursor_ion_mass": float(spectra[5]), "peaks": [],
                                                     "accession": hmdb_id, "exact_mass": spectra[2],
                                                     "actual_accession": spectra[4]}

            data_categories[name_msp][hmdb_id]["exact_mass"] = float(spectra[2])

            data_categories[name_msp][hmdb_id]["mol"] = mol
            data_categories[name_msp][hmdb_id]["smiles"] = Chem.MolToSmiles(mol)

            data_categories[name_msp][hmdb_id]["chemical_formula"] = get_elements(mol)

            chemical_formula = []
            for element in ["C", "H", "N", "O", "P", "S"]:
                data_categories[name_msp][hmdb_id][element] = data_categories[name_msp][hmdb_id]["chemical_formula"][element]
                chemical_formula.append(data_categories[name_msp][hmdb_id]["chemical_formula"][element])

            data_categories[name_msp][hmdb_id]["mc"] = chemical_formula
            data_categories[name_msp][hmdb_id]["chemical_formula"] = ""

            data_categories[name_msp][hmdb_id]["peaks"] = spectra[6]
            data_categories[name_msp][hmdb_id]["intensities"] = spectra[7]

        msp_db.close()

    return data_categories
