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

import numpy
from math import sqrt


def find_path(mass_list, sum_matrix, n, mass, max_subset_length, path=None):
    """
    Recursive solution for backtracking through the dynamic programming boolean matrix. All possible subsets are found

    :param mass_list: A list of masses from which to identify subsets.

    :param mass: The target mass of the sum of the substructures.

    :param sum_matrix: The dynamic programming boolean matrix.

    :param n: The size of mass_list.

    :param max_subset_length: The maximum length of subsets to return. Allows the recursive backtracking algorithm to
        terminate early in many cases, significantly improving runtime.

    :param path: List for keeping track of the current subset.

    :return: Generates of lists containing the masses of valid subsets.
    """

    if path is None:
        path = []

    # base case - the path has generated a correct solution
    if mass == 0:
        yield sorted(path)
        return

    # stop running when we overshoot the mass
    elif mass < 0:
        return

    # can we sum up to the target value using the remaining masses? recursive call
    elif sum_matrix[n][mass]:
        yield from find_path(mass_list, sum_matrix, n - 1, mass, max_subset_length, path)

        if len(path) < max_subset_length:
            path.append(mass_list[n-1])

            yield from find_path(mass_list, sum_matrix, n - 1, mass - mass_list[n - 1], max_subset_length, path)
            path.pop()


def subset_sum(mass_list, mass, max_subset_length=3):
    """
    Dynamic programming implementation of subset sum. Note that, whilst this algorithm is pseudo-polynomial, the
    backtracking algorithm for obtaining all possible subsets has exponential complexity and so remains unsuitable
    for large input values.  This does, however, tend to perform a lot better than non-sum_matrix implementations, as
    we're no longer doing sums multiple times and we've cut down the operations performed during the exponential portion
    of the method.

    :param mass_list: A list of masses from which to identify subsets.

    :param mass: The target mass of the sum of the substructures.

    :param max_subset_length: The maximum length of subsets to return. Allows the recursive backtracking algorithm to
        terminate early in many cases, significantly improving runtime.

    :return: Generates of lists containing the masses of valid subsets.
    """

    n = len(mass_list)

    # initialise dynamic programming array
    sum_matrix = numpy.ndarray([n + 1, mass + 1], bool)

    # subsets can always equal 0
    for i in range(n+1):
        sum_matrix[i][0] = True

    # empty subsets do not have non-zero sums
    for i in range(mass):
        sum_matrix[0][i + 1] = False

    # fill in the remaining boolean matrix
    for i in range(n):
        for j in range(mass+1):
            if j >= mass_list[i]:
                sum_matrix[i + 1][j] = sum_matrix[i][j] or sum_matrix[i][j - mass_list[i]]
            else:
                sum_matrix[i + 1][j] = sum_matrix[i][j]

    # backtrack through the matrix recursively to obtain all solutions
    return find_path(mass_list, sum_matrix, n, mass, max_subset_length)


def cosine_spectrum_similarity(real_mzs, candidate_mzs):
    """
    Database fragmentation scoring based on the cosine similarity method. Adapted for the lack of intensities
    available for the candidate compound.

    :param real_mzs: The mz values for the original MSn spectrum.

    :param candidate_mzs: The theoretical mz values for a candidate compound. Should have the same order as `real_mzs`
        and should have a value of `0` when there is no match for the candidate for a peak in the original spectrum.

    :return: Similarity metric for the two spectra.
    """

    # get weighted peaks
    real_weighted = [(mz ** 2) for mz in real_mzs]
    candidate_weighted = [(mz ** 2) for mz in candidate_mzs]

    def dot(E, D):
        return sum(e * d for e, d in zip(E, D))

    def cosine_similarity(E, D):
        return dot(E, D) / (sqrt(dot(E, E)) * sqrt(dot(D, D)))

    return cosine_similarity(real_weighted, candidate_weighted)
