#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019-2020 Ralf Weber
#
# This file is part of MetaboBlend.
#
# MetaboBlend is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MetaboBlend is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MetaboBlend.  If not, see <https://www.gnu.org/licenses/>.
#

import os
import csv
import sqlite3
from math import sqrt


class ResultsDb:
    """
    Methods for interacting with the SQLITE3 results database, as created by
    :py:meth:`metaboblend.build_structures.annotate_msn`.

    :param path_results: Directory to which results will be written.
    """

    def __init__(self, path_results, msn=True):
        """Constructor method."""

        self.path_results = path_results
        self.path_results_db = os.path.join(self.path_results, "metaboblend_results.sqlite")
        self.msn = msn

        self.conn = None
        self.cursor = None
        self.open()

        self.substructure_combo_id = 0

    def open(self):
        """ Opens connection to the SQLite3 database and creates custom functions. """

        self.conn = sqlite3.connect(self.path_results_db)
        self.cursor = self.conn.cursor()

        # MS-FINDER method of calculating bde scores
        self.conn.create_function("CALC_BDE_SCORE", 2, lambda bde, max_bde: sqrt(1 - (bde / max_bde)))

        def calc_result_score(even, bde_score):

            # doesn't follow hydrogenation rules (MS-FINDER)
            if even == 0:
                bde_score *= 0.5

            # add some base weight to peaks such that: scores >= 0.5
            bde_score = (0.5 + bde_score) / 1.5

            return bde_score

        self.conn.create_function("CALC_RESULT_SCORE", 2, calc_result_score)

    def create_results_db(self):
        """ Generates a new results database. """

        self.close()

        if os.path.exists(self.path_results_db):
            os.remove(self.path_results_db)

        self.open()

        self.cursor.execute("""CREATE TABLE queries (
                                   ms_id_num INTEGER PRIMARY KEY,
                                   ms_id TEXT,
                                   exact_mass NUMERIC,
                                   C INTEGER,
                                   H INTEGER,
                                   N INTEGER,
                                   O INTEGER,
                                   P INTEGER,
                                   S INTEGER,
                                   ppm INTEGER,
                                   ha_min INTEGER,
                                   ha_max INTEGER,
                                   max_atoms_available INTEGER,
                                   max_degree INTEGER,
                                   max_n_substructures INTEGER,
                                   hydrogenation_allowance INTEGER,
                                   isomeric_smiles INTEGER)""")

        if self.msn:
            self.cursor.execute("""CREATE TABLE spectra (
                                       ms_id_num INTEGER,
                                       fragment_id INTEGER,
                                       neutral_mass NUMERIC,
                                       max_bde NUMERIC,
                                       PRIMARY KEY (ms_id_num, fragment_id))""")

        self.create_structures_table()

        self.cursor.execute("""CREATE TABLE results (
                                   ms_id_num INTEGER,
                                   fragment_id INTEGER,
                                   structure_smiles TEXT,
                                   even BOOLEAN,
                                   result_score NUMERIC,
                                   PRIMARY KEY(ms_id_num, fragment_id, structure_smiles))""")

        self.cursor.execute("""CREATE TABLE substructure_combos (
                                           substructure_combo_id INTEGER,
                                           ms_id_num INTEGER,
                                           fragment_id INTEGER,
                                           structure_smiles TEXT,
                                           bde INTEGER,
                                           bde_score NUMERIC,
                                           PRIMARY KEY (substructure_combo_id))""")

        self.cursor.execute("""CREATE TABLE substructures (
                                           substructure_combo_id INTEGER,
                                           substructure_position_id INTEGER,
                                           substructure_smiles TEXT,
                                           PRIMARY KEY (substructure_combo_id, substructure_position_id))""")

        self.conn.commit()

    def create_structures_table(self):
        """ Create structures table. """

        self.cursor.execute("""CREATE TABLE structures (
                                   ms_id_num INTEGER,
                                   structure_smiles TEXT,
                                   frequency INTEGER,
                                   frequency_score NUMERIC,
                                   PRIMARY KEY (ms_id_num, structure_smiles))""")

    def add_ms(self, msn_data, ms_id, ms_id_num, parameters):
        """
        Add entries to the `queries` and `spectra` tables.

        :param msn_data: Dictionary in the form
            `msn_data[id] = {mf: [C, H, N, O, P, S], exact_mass: float, fragment_masses: []}`. id represents a unique
            identifier for a given spectral tree or fragmentation spectrum, mf is a list of integers referring to the
            molecular formula of the structure of interest, exact_mass is the mass of this molecular formula to >=4d.p.
            and fragment_masses are neutral fragment masses generated by this structure used to inform candidate
            scoring. See :py:meth:`metaboblend.build_structures.annotate_msn`.

        :param ms_id: Unique identifier for the annotation of a single metabolite.

        :param ms_id_num: Unique numeric identifier for the annotation of a single metaoblite.

        :param parameters: List of parameters, in the form: [ppm, ha_min, ha_max, max_atoms_available, max_degree,
            max_n_substructures, hydrogenation_allowance, isomeric_smiles]. See
            :py:meth:`metaboblend.build_structures.annotate_msn`.
        """

        for i, parameter in enumerate(parameters):
            if parameter is None:
                parameters[i] = "NULL"
            elif isinstance(parameter, bool):
                parameters[i] = int(parameter)

        self.cursor.execute("""INSERT INTO queries (
                                   ms_id,
                                   ms_id_num,
                                   exact_mass,
                                   C, H, N, O, P, S,
                                   ppm,
                                   ha_min,
                                   ha_max,
                                   max_atoms_available,
                                   max_degree,
                                   max_n_substructures,
                                   hydrogenation_allowance,
                                   isomeric_smiles
                               ) VALUES ('{}', {}, {}, '{}', '{}', '{}', '{}', '{}', '{}', {})""".format(
                                   ms_id,
                                   ms_id_num,
                                   msn_data[ms_id]["exact_mass"],
                                   msn_data[ms_id]["mf"][0], msn_data[ms_id]["mf"][1],
                                   msn_data[ms_id]["mf"][2], msn_data[ms_id]["mf"][3],
                                   msn_data[ms_id]["mf"][4], msn_data[ms_id]["mf"][5],
                                   ", ".join([str(p) for p in parameters])
                               ))

        self.conn.commit()

    def add_results(self, ms_id_num, smi_dict, fragment_mass=None, fragment_id=None, retain_substructures=False):
        """
        Record which smiles were generated for a given fragment mass.

        :param ms_id_num: Unique identifier for the annotation of a single metabolite.

        :param smi_dict: The fragment and substructure smiles generated by the annotation of a single peak for a single
            metabolite.

        :param fragment_mass: The neutral fragment mass that has been annotated.

        :param fragment_id: The unique identifier for the fragment mass that has been annotated.

        :param retain_substructures: If True, record substructures in the results DB.
        """

        # if annotating msn spectra, fill the spectra table
        if self.msn:

            # get the maximum BDE across all structures generated for this particular fragment ion
            max_bde = 0

            for structure_smiles in smi_dict.keys():
                max_bde = max(max_bde, max(smi_dict[structure_smiles]["bdes"]))

            self.cursor.execute("""INSERT OR IGNORE INTO spectra (
                                       ms_id_num,
                                       fragment_id,
                                       neutral_mass,
                                       max_bde
                                   ) VALUES ('{}', {}, {}, {})
                                """.format(
                                       ms_id_num,
                                       fragment_id,
                                       fragment_mass,
                                       max_bde
                                   ))
        else:
            fragment_id = "NULL"

        # unique smiles candidates for this fragment
        for structure_smiles in smi_dict.keys():

            if not self.msn:
                smi_dict[structure_smiles]["even"] = 1

            # for each combination of substructures that generated the candidate
            for i in range(len(smi_dict[structure_smiles]["substructures"])):

                self.cursor.execute("""INSERT INTO substructure_combos (
                                           substructure_combo_id,
                                           ms_id_num,
                                           fragment_id,
                                           structure_smiles,
                                           bde
                                       ) VALUES ({}, {}, {}, '{}', {})
                                    """.format(
                                           self.substructure_combo_id,
                                           ms_id_num,
                                           fragment_id,
                                           structure_smiles,
                                           smi_dict[structure_smiles]["bdes"][i]
                                    ))

                if retain_substructures:
                    for j, substructure in enumerate(smi_dict[structure_smiles]["substructures"][i]):

                        self.cursor.execute("""INSERT INTO substructures (
                                                   substructure_combo_id,
                                                   substructure_position_id,
                                                   substructure_smiles
                                               ) VALUES ({}, {}, '{}')
                                            """.format(
                                                   self.substructure_combo_id,
                                                   j,
                                                   substructure
                                            ))

                self.substructure_combo_id += 1

            self.cursor.execute("""INSERT INTO results (
                                       ms_id_num,
                                       fragment_id,
                                       structure_smiles,
                                       even
                                   ) VALUES ({}, {}, '{}', {})
                                """.format(
                                       ms_id_num,
                                       fragment_id,
                                       structure_smiles,
                                       int(smi_dict[structure_smiles]["even"])
                                ))

        self.conn.commit()

    def calculate_scores(self, ms_id_num):
        """
        Scores cannot be calculated while generating the various tables. Must be completed after the entire structure
        generation process of a metabolite. For instance, the maximum BDE across all sets of substructures generated
        for a metabolite can only be ascertained once all these sets of substructures have been recorded.

        Does the calculations for aggregating structure candidate scores within SQL by updating columns that have been
        ignored thus far. More complex calculations are written in python as SQL functions, as defined in
        py:meth:`metaboblend.results.open`.

        :param ms_id_num: Unique identifier for the annotation of a single metabolite.
        """

        if not self.msn:
            self.cursor.execute("""INSERT INTO structures (ms_id_num, structure_smiles, frequency)
                                       SELECT ms_id_num, structure_smiles, COUNT(*)
                                           FROM results
                                           WHERE ms_id_num = {}
                                           GROUP BY structure_smiles""".format(ms_id_num))

            return

        self.cursor.execute("SELECT COUNT(*), MAX(max_bde) FROM spectra WHERE ms_id_num = %s" % ms_id_num)
        num_fragments, max_bde = list(self.cursor.fetchall())[0]

        # calculate the BDE score for each combination of substructures
        self.cursor.execute("""UPDATE substructure_combos
                                   SET bde_score = CALC_BDE_SCORE(bde, {})
                                   WHERE ms_id_num = {}
                            """.format(max_bde, ms_id_num))

        # aggregate bde scores for each peak/candidate structure
        # updates results by aggregating the BDE scores (selects the max score for the peak/candidate structure)
        # then correlating the result of this query with the results table
        self.cursor.execute("""UPDATE results
                                   SET result_score = CALC_RESULT_SCORE(even, (
                                       SELECT max_bde_score FROM (                                   
                                           SELECT fragment_id, structure_smiles, MAX(bde_score) AS max_bde_score
                                               FROM substructure_combos
                                               WHERE ms_id_num = {} 
                                               GROUP BY fragment_id, structure_smiles
                                           )
                                           WHERE fragment_id = results.fragment_id 
                                               AND structure_smiles = results.structure_smiles
                                       ))
                                   WHERE ms_id_num = {}
                            """.format(ms_id_num, ms_id_num))

        # aggregate results scores across the spectrum for each unique structure candidate
        self.cursor.execute("""INSERT INTO structures (ms_id_num, structure_smiles, frequency, frequency_score)
                                   SELECT ms_id_num, structure_smiles, COUNT(*), SUM(result_score) / {}
                                       FROM results
                                       WHERE ms_id_num = {}
                                       GROUP BY structure_smiles""".format(num_fragments, ms_id_num))

        self.conn.commit()

    def recalculate_scores(self):
        """ Re-calculates scores for the results DB. """

        self.cursor.execute("DROP TABLE IF EXISTS structures")
        self.create_structures_table()

        self.cursor.execute("SELECT DISTINCT ms_id_num FROM queries")
        ms_id_nums = [row[0] for row in self.cursor.fetchall()]

        for i in ms_id_nums:
            self.calculate_scores(i)

    def get_structures(self, ms_id_num):
        """
        Gets smiles of generated structures. In the case of the MSn annotation workflow, also gets structure
        frequencies.

        :param ms_id_num: Unique identifier for the annotation of a single metabolite.

        :return: In the case of simple structure generation, returns a set of smiles strings for output structures.
            For the MSn annotation workflow, returns a dictionary with smiles as keys and the number of peaks for which
            the smiles were generated as values.
        """

        if self.msn:
            msn_str = ", frequency"
        else:
            msn_str = ""

        self.cursor.execute("""SELECT structure_smiles{} FROM structures 
                                   WHERE ms_id_num = {}
                            """.format(msn_str, ms_id_num))

        if self.msn:
            return [t for t in self.cursor.fetchall()]
        else:
            return [item for t in self.cursor.fetchall() for item in t]

    def generate_csv_output(self):
        """ Generate CSV file output for i) queries and tool parameters and ii) structures generated. """

        with open(os.path.join(self.path_results, "metaboblend_queries.csv"), "w", newline="") as results_file, \
                open(os.path.join(self.path_results, "metaboblend_structures.csv"), "w", newline="") as ms_file:

            results_writer = csv.writer(results_file, delimiter=",")
            ms_writer = csv.writer(ms_file, delimiter=",")

            results_writer.writerow(["ms_id", "exact_mass", "C", "H", "N", "O", "P", "S", "ppm", "ha_min", "ha_max",
                                     "max_atoms_available", "max_degree", "max_n_substructures",
                                     "hydrogenation_allowance", "isomeric_smiles"])

            self.cursor.execute("SELECT * FROM queries")

            for query in self.cursor.fetchall():
                results_writer.writerow(query)

            ms_writer.writerow(["ms_id", "smiles", "frequency", "exact_mass", "C", "H", "N", "O", "P", "S"])

            self.cursor.execute("SELECT * FROM structures")

            for structure in self.cursor.fetchall():
                ms_writer.writerow(structure)

    def close(self):
        """ Close the connection to the SQLITE3 database. """

        self.conn.close()
