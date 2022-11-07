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

import argparse
import os

from metaboblend import __version__

from . import databases
from . import build_structures


def main():  # pragma: no cover

    print(("Executing metaboblend version %s." % __version__))

    parser = argparse.ArgumentParser(description='Python package for de novo structural elucidation of small molecules in mass spectrometry-based Metabolomics',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers(dest='step')

    parser_a = subparsers.add_parser('to-add', help='to-add')

    parser_a.add_argument('-i', '--input',
                          type=str, required=True,
                          help="to-add")

    args = parser.parse_args()

    print(args)

    if args.step == "to-add":
        pass


if __name__ == "__main__":
    main()
