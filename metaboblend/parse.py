#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019-2020 Jack Gisby, Ralf Weber
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

import io
import re
import sys
import copy
import warnings
from collections import OrderedDict
import xml.etree.ElementTree as ElementTree


def parse_ms_data(ms_data, msn=True):
    """
    Parse raw data provided by user and yield formatted input data. Decides what type of data has been provided
    (i.e. whether a dictionary has been given vs path to MSP file; if a dictionary, checks whether neutral masses
    need to be calculated from precursor ions).

    :param ms_data: Dictionary containing input data or path to an MSP file.

    :param msn: If True, formats the data for use by :py:meth:`metaboblend.build_structures.annotate_msn`; else, formats
        input data for use by :py:meth:`metaboblend.build_structures.generate_structures`. Only relevant if a
        dictionary has been provided.

    :return: Yields a dictionary for use by build functions to generate structures.
    """

    if isinstance(ms_data, dict):
        for i, ms_id in enumerate(ms_data.keys()):

            ms_data[ms_id]["ms_id"] = ms_id

            # check if user has provided a neutralised mass or ionised mz values
            if "neutral_fragment_masses" in ms_data[ms_id].keys() and "exact_mass" in ms_data[ms_id].keys():
                which = "none"

            elif "exact_mass" in ms_data[ms_id].keys():
                if msn:
                    which = "fragments"
                else:
                    which = "none"

            elif "neutral_fragment_masses" in ms_data[ms_id].keys() or not msn:
                which = "precursor"

            else:
                which = "both"

            yield precursor_ions_to_neutral_masses(ms_data[ms_id], which)

    else:
        yield from parse_msp(ms_data)


def precursor_ion_to_neutral_mass(mass, precursor_type, get_mode=False):
    """
    Convert precursor ion to predicted neutral mass for substructure searching.

    :param mass: Charged mass to be neutralised.

    :param precursor_type: Type of precursor ion.

    :param get_mode: If True, also return the ion mode (positive or negative).

    :return: Neutral mass.
    """

    # conversions
    precursor_dict = {"[M+H]+":  1.007276,
                      "[M+Na]+": 22.989221,
                      "[M+K]+": 38.963158,
                      "[M-H]-": -1.007276,
                      "[M+Cl]-": 34.969401,
                      "[M+Na-2H]-": 20.974668,
                      "[M+K-2H]-": 36.948605,
                      "[M+Hac-H]-": 59.013853}

    precursor_mode = {"[M+H]+":  "+", "[M+Na]+": "+", "[M+K]+": "+",
                      "[M-H]-": "-", "[M+Cl]-": "-", "[M+Na-2H]-": "-", "[M+K-2H]-": "-", "[M+Hac-H]-": "-"}

    if get_mode:
        return mass - precursor_dict[precursor_type], precursor_mode[precursor_type]
    else:
        return mass - precursor_dict[precursor_type]


def precursor_ions_to_neutral_masses(ms_dict, which="both"):
    """
    Convert precursor ion and fragment ions to neutral.

    :param ms_dict: Dictionary used by build functions to generate structures. Converts the precursor ion mass and/or
        the fragment ions to their respective neutral masses.

    :param which: Whether to convert the precursor ion ("precursor"), the fragment ions ("fragments") or both ("both")
        to their respective neutral masses. If which is "none", returns the original dictionary.

    :return: Returns `ms_dict` with additional items corresponding to neutralised masses.
    """

    if which == "precursor" or which == "both":
        ms_dict["exact_mass"], ms_dict["ion_mode"] = precursor_ion_to_neutral_mass(ms_dict["precursor_mz"],
                                                                                   ms_dict["precursor_type"],
                                                                                   get_mode=True)

    if which == "fragments" or which == "both":

        ms_dict["neutral_fragment_masses"] = []

        for fragment_ion_mass in ms_dict["fragment_mzs"]:
            ms_dict["neutral_fragment_masses"].append(precursor_ion_to_neutral_mass(fragment_ion_mass,
                                                                                    ms_dict["precursor_type"]))

    return ms_dict


def parse_msp(msp_path):
    """
    Parse msp files and yield data for each compound. Accepts MSP files in MoNa or MassBank format. We expect that
    the following are provided in the MSP:

    - A unique accession ID.
    - The molecular formula of the compound.
    - The precursor mz representing the mass of the charged precursor ion.
    - Fragment mzs representing masses of charged fragment ions.
    - The type of precursor, e.g. "[M+H]+".

    Code adapted from `msp2db` (https://github.com/computational-metabolomics/msp2db/blob/master/msp2db/parse.py).

    :param msp_path: Path of an MSP file to be converted into a dictionary.

    :return: Dictionary in a form useable by :py:meth:`metaboblend.build_structures.annotate_msn` and
        :py:meth:`metaboblend.build_structures.generate_structures`.
    """

    meta_parse = get_msp_regex()
    reached_spectra = False

    empty_dict = {"ms_id": None, "mf": None, "precursor_mz": None, "precursor_type": None, "fragment_mzs": []}
    entry_dict = copy.deepcopy(empty_dict)

    with open(msp_path, "r") as msp_file:

        for line in msp_file:

            line = re.sub('^(.{2}\\$)', "", line)  # remove "XX$" from line start in massbank files

            if reached_spectra:
                if line in ["\n", "\r\n", "//\n", "//\r\n", "", "//"]:  # reached end of spectra

                    yield reformat_msp_input(entry_dict)  # completed entry ready for sending to build

                    entry_dict = copy.deepcopy(empty_dict)
                    reached_spectra = False

                else:  # add peak
                    entry_dict["fragment_mzs"].append(float(line.split()[0]))

            else:
                for meta_type in meta_parse.keys():
                    for meta_re in meta_parse[meta_type]:

                        re_query = re.search(meta_re, line, re.IGNORECASE)

                        if re_query:
                            entry_dict[meta_type] = re_query.group(1).strip()

                if re.match("^Num Peaks(.*)$", line, re.IGNORECASE) or re.match("^PEAK:(.*)", line, re.IGNORECASE):
                    reached_spectra = True  # reached line prior to spectra

    if entry_dict != empty_dict:
        yield reformat_msp_input(entry_dict)


def reformat_msp_input(entry_dict):
    """
    Reformat input for use by build functions.

    :param entry_dict: Dictionary containing MSn information extracted from an MSP file (by
        :py:meth:`metaboblend.parse.parse_msp`. The dictionary must contain the following:

        - ms_id - a unique accession number
        - mf - the molecular formula of the compound (in the format "CXHXNXOXPXSX")
        - precursor_mz - mz representing the mass of the charged precursor ion
        - precursor_type - the type of precursor ion (e.g. "[M+H]+")
        - fragment_mzs - mz(s) representing the mass of charged fragment ions

    :return: If the correct inputs were not provided in the MSP (and, hence, were not available in `entry_dict`),
        returns None (and generates a warning with i) the accession (if available) and ii) the variable that was not
        able to be extracted from the MSP). Else, returns the same dictionary after reformatting the molecular formula,
        using :py:meth:`metaboblend.parse.mc_to_list`, and converting the precursor ions to their corresponding
        neutral masses.
    """

    if entry_dict["mf"] is not None:  # convert from C5H6... to [5, 6, ...]
        entry_dict["mf"] = mc_to_list(entry_dict["mf"])

    for key in ["ms_id", "mf", "precursor_mz", "precursor_type"]:  # required for the tool to function
        if entry_dict[key] is None:
            if key == "ms_id":
                warnings.warn("Entry ignored from MSP file due to lack of accession in MSP file")
            else:
                warnings.warn("Entry " + entry_dict["ms_id"] + " removed due to lack of valid " + key + " in MSP file")
            return None

    entry_dict["precursor_mz"] = float(entry_dict["precursor_mz"])

    if len(entry_dict["fragment_mzs"]) == 0:  # require a spectra to annotate
        warnings.warn("No fragments were identified for " + entry_dict["ms_id"] + " in MSP file")
        return None

    return precursor_ions_to_neutral_masses(entry_dict)


def mc_to_list(mc):
    """
    Convert molecular formula string to list format.

    :param mc: Molecular formula (in the format "C1H2N3O4P5S6")

    :return: Molecular formula (in the format `[1, 2, 3, 4, 5, 6]`)
    """

    if isinstance(mc, list):
        return mc

    mc_list = [0, 0, 0, 0, 0, 0]
    element_positions = {"C": 0, "H": 1, "N": 2, "O": 3, "P": 4, "S": 5}

    # seperates out the formula into [letter, number, letter, number, ...]
    mc = re.findall(r"[A-Z][a-z]*|\d+", re.sub("[A-Z][a-z]*(?![\da-z])", r"\g<0>1", mc))

    for i, substring in enumerate(mc):

        if i % 2 == 0:  # in case of letter
            try:
                element_position = element_positions[substring]
            except KeyError:  # element not in C, H, N, O, P, S
                return None

        else:  # record number following the letter
            mc_list[element_position] = int(substring)

    return mc_list


def get_msp_regex():
    """ Dictionary of regular expressions for parsing msp metadata. """

    meta_parse = {"ms_id":          ["^accession(?:=|:)(.*)$", "^DB#(?:=|:)(.*)$", "^ACCESSION:(.*)$"],  # use accession as ms_id
                  "mf":             ["^molecular formula(?:=|:)(.*)$", "^formula:(.*)$"],
                  "precursor_type": ["^precursor.*type(?:=|:)(.*)$", "^adduct(?:=|:)(.*)$", "^MS\$FOCUSED_ION:\s+PRECURSOR_TYPE\s+(.*)$"],
                  "precursor_mz":   ["^precursor m/z(?:=|:)\s*(\d*[.,]?\d*)$", "^precursor.*mz(?:=|:)\s*(\d*[.,]?\d*)$", "^MS\$FOCUSED_ION:\s+PRECURSOR_M/Z\s+(\d*[.,]?\d*)$"]}

    return meta_parse


def reformat_xml(source, encoding="utf8"):
    """
    Reformats HMDB xml files to be compatible with :py:meth:`metaboblend.databases.parse_xml`; some such files do not
    contain a `<hmdb xmlns="http://www.hmdb.ca">` header.

    :param source: Path to file to be reformatted.

    :param encoding: Encoding of source file.

    :return: Source file destination.
    """

    with io.open(source, "r", encoding=encoding) as xml:
        xml_contents = xml.readlines()
        if "hmdb" in xml_contents[1]:
            return source

        xml_contents.insert(1, "<hmdb xmlns=\"http://www.hmdb.ca\"> \n")

    with io.open(source, "w", encoding=encoding) as xml:
        xml_contents = "".join(xml_contents)
        xml.write(xml_contents)
        xml.write("</hmdb>")

    return source


def parse_xml(source, encoding="utf8", reformat=False):
    """
    Parses the contents of HMDB xml files to to extract information for the generation of substructures.

    :param source: Source file destination.

    :param encoding: Encoding of source file.

    :param reformat: Whether to apply :py:meth:`metaboblend.databases.reformat_xml` to the XML file. Is required for
        XML files recording single metabolites.

        * **True** Add a `<hmdb xmlns="http://www.hmdb.ca">` header to the XML file before parsing.

        * **False** Parse the XML file as it is (recommended if header is present).

    :return: The XML file converted to a dictionary.
    """

    if reformat:
        reformat_xml(source, encoding)

    with io.open(source, "r", encoding=encoding) as inp:
        record_out = OrderedDict()

        inp.readline()
        inp.readline()

        xml_record = ""
        path = []

        for line in inp:
            xml_record += line
            if line == "</metabolite>\n" or line == "</drug>\n":

                if sys.version_info[0] == 3:
                    inp = io.StringIO(xml_record)
                else:
                    inp = io.BytesIO(xml_record.encode('utf-8').strip())

                for event, elem in ElementTree.iterparse(inp, events=("start", "end")):
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
