#!/usr/bin/env python3
"""Exact finite checks of cyclic second-family states and shield identities."""

from collections import Counter
from fractions import Fraction
from itertools import product
import json
from math import prod


class Quadratic:
    def __init__(self, dimension, rational=0, radical=0):
        self.dimension = dimension
        self.rational = Fraction(rational)
        self.radical = Fraction(radical)

    def __add__(self, other):
        if not isinstance(other, Quadratic):
            other = Quadratic(self.dimension, other)
        assert self.dimension == other.dimension
        return Quadratic(self.dimension, self.rational + other.rational,
                         self.radical + other.radical)

    def __mul__(self, other):
        if not isinstance(other, Quadratic):
            other = Quadratic(self.dimension, other)
        assert self.dimension == other.dimension
        return Quadratic(self.dimension,
                         self.rational * other.rational
                         + self.dimension * self.radical * other.radical,
                         self.rational * other.radical + self.radical * other.rational)

    __radd__ = __add__
    __rmul__ = __mul__

    def __neg__(self):
        return Quadratic(self.dimension, -self.rational, -self.radical)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    def __eq__(self, other):
        if not isinstance(other, Quadratic):
            other = Quadratic(self.dimension, other)
        return (self.dimension, self.rational, self.radical) == (
            other.dimension, other.rational, other.radical)

    def __bool__(self):
        return bool(self.rational or self.radical)


def add_entry(matrix, row, column, value):
    position = row, column
    updated = matrix.get(position, 0) + value
    if updated:
        matrix[position] = updated
    else:
        matrix.pop(position, None)


def transpose_bob(matrix):
    output = {}
    for (row, column), value in matrix.items():
        row_a, row_b, row_shield_a, row_shield_b = row
        column_a, column_b, column_shield_a, column_shield_b = column
        new_row = row_a, column_b, row_shield_a, column_shield_b
        new_column = column_a, row_b, column_shield_a, row_shield_b
        add_entry(output, new_row, new_column, value)
    return output


def add_projector(matrix, support, coefficient):
    for row, column in product(support, repeat=2):
        add_entry(matrix, row, column, coefficient)


def check_group(moduli):
    labels = list(product(*(range(modulus) for modulus in moduli)))
    dimension = prod(moduli)

    def add(left, right):
        return tuple((first + second) % modulus
                     for first, second, modulus in zip(left, right, moduli))

    def subtract(left, right):
        return tuple((first - second) % modulus
                     for first, second, modulus in zip(left, right, moduli))

    zero = tuple(0 for _ in moduli)
    for shift in labels:
        assert {add(row, shift) for row in labels} == set(labels)
    multiplicities = Counter(subtract(column, row) for row, column in product(labels, repeat=2))
    assert set(multiplicities.values()) == {dimension}
    entries_checked = 0
    wrong_same_shield_entries = 0
    for row_a, row_b, column_a, column_b in product(labels, repeat=4):
        partial_transpose_sum = int(add(row_a, column_b) == add(column_a, row_b))
        difference_projector = int(subtract(row_b, row_a) == subtract(column_b, column_a))
        assert partial_transpose_sum == difference_projector
        same_shield_transpose = int(subtract(column_b, row_a) == subtract(row_b, column_a))
        wrong_same_shield_entries += same_shield_transpose != difference_projector
        entries_checked += 1
    exponent_two = all(add(label, label) == zero for label in labels)
    assert (wrong_same_shield_entries == 0) == exponent_two
    return {'moduli': moduli, 'dimension': dimension,
            'shield_entries_checked': entries_checked,
            'sum_to_difference_identity': True,
            'same_shield_identity': exponent_two,
            'same_shield_mismatches': wrong_same_shield_entries}


def check_state(dimension):
    labels = range(dimension)
    root = Quadratic(dimension, 0, 1)
    weight_p = Quadratic(dimension, Fraction(dimension, dimension - 1),
                         Fraction(-1, dimension - 1))
    weight_q = Quadratic(dimension, Fraction(-1, dimension - 1),
                         Fraction(1, dimension - 1))
    assert weight_p + weight_q == 1
    assert root * weight_q == weight_p
    state = {}
    for row, column in product(labels, repeat=2):
        support = [(0, 0, row, column), (1, 1, column, (row + column) % dimension)]
        add_projector(state, support, weight_p * Fraction(1, 2 * dimension ** 2))
    for sector in labels:
        support = [(0, 1, row, (sector - row) % dimension) for row in labels]
        add_projector(state, support, weight_q * Fraction(1, 2 * dimension ** 2))
        add_entry(state, (1, 0, sector, sector), (1, 0, sector, sector),
                  weight_q * Fraction(1, 2 * dimension))
    predicted_transpose = {}
    for row, column in product(labels, repeat=2):
        for key in [(0, 0), (1, 1)]:
            basis = key + (row, column)
            add_entry(predicted_transpose, basis, basis,
                      weight_p * Fraction(1, 2 * dimension ** 2))
    for sector in labels:
        support = [(0, 1, row, (row + sector) % dimension) for row in labels]
        add_projector(predicted_transpose, support,
                      weight_q * Fraction(1, 2 * dimension ** 2))
        diagonal = (1, 0, sector, sector)
        add_entry(predicted_transpose, diagonal, diagonal,
                  weight_q * Fraction(1, 2 * dimension))
        for basis in support:
            coefficient = weight_q * root * Fraction(1, 2 * dimension ** 2)
            add_entry(predicted_transpose, basis, diagonal, coefficient)
            add_entry(predicted_transpose, diagonal, basis, coefficient)
    actual_transpose = transpose_bob(state)
    assert actual_transpose == predicted_transpose
    assert transpose_bob(actual_transpose) == state
    assert sum(value for (row, column), value in state.items() if row == column) == 1
    assert sum(value for (row, column), value in actual_transpose.items() if row == column) == 1
    assert all(state.get((column, row), 0) == value
               for (row, column), value in state.items())
    occupied_z = [(column, (row + column) % dimension)
                  for row, column in product(labels, repeat=2)]
    assert len(set(occupied_z)) == dimension ** 2
    sum_supports = [set((row, (sector - row) % dimension) for row in labels)
                    for sector in labels]
    difference_supports = [set((row, (row + sector) % dimension) for row in labels)
                           for sector in labels]
    for first, second in product(labels, repeat=2):
        expected_overlap = dimension if first == second else 0
        assert len(sum_supports[first] & sum_supports[second]) == expected_overlap
        assert len(difference_supports[first] & difference_supports[second]) == expected_overlap
    mean_h = sum(value * (2 - 2 * row[0] - 2 * row[1])
                 for (row, column), value in state.items() if row == column)
    mean_h_squared = sum(value * (2 - 2 * row[0] - 2 * row[1]) ** 2
                         for (row, column), value in state.items() if row == column)
    assert mean_h == 0
    assert mean_h_squared == weight_p * 4
    fisher_information = sum((weight_p * Fraction(1, dimension ** 2)) * 16
                             for _ in range(dimension ** 2))
    assert fisher_information == weight_p * 16
    noise_checks = 0
    for visibility in [Fraction(0), Fraction(1, 5), Fraction(1, 2), Fraction(1)]:
        noise_floor = Fraction(1 - visibility, 4 * dimension ** 2)
        eigenvalue_plus = weight_p * Fraction(visibility, dimension ** 2) + noise_floor
        eigenvalue_minus = Quadratic(dimension, noise_floor)
        difference = eigenvalue_plus - eigenvalue_minus
        denominator = eigenvalue_plus + eigenvalue_minus
        numerator = difference * difference * (16 * dimension ** 2)
        alternative_numerator = weight_p * weight_p * (32 * visibility ** 2)
        alternative_denominator = weight_p * (2 * visibility) + 1 - visibility
        assert numerator * alternative_denominator == alternative_numerator * denominator
        noise_checks += 1
    return {'dimension': dimension, 'state_nonzero_entries': len(state),
            'partial_transpose_nonzero_entries': len(actual_transpose),
            'density_and_transpose_identity': True,
            'state_rank': dimension ** 2 + 2 * dimension,
            'partial_transpose_rank': 2 * dimension ** 2 + dimension,
            'orthogonality_and_normalization': True,
            'fisher_information': '16*sqrt(d)/(1+sqrt(d))',
            'noise_parameter_checks': noise_checks}


def main():
    cyclic = [check_group([dimension]) for dimension in range(2, 17)]
    groups = [check_group(moduli) for moduli in [[2, 2], [2, 3], [2, 2, 2], [3, 3], [2, 2, 2, 2]]]
    states = [check_state(dimension) for dimension in range(2, 17)]
    assert cyclic[1]['same_shield_mismatches'] > 0
    print(json.dumps({
        'arithmetic': 'exact rational coefficients in Q[r]/(r^2-d)',
        'coverage': 'Finite cyclic dimensions 2 through 16 and five product groups; the universal statements have written proofs in note.tex.',
        'status': 'pass', 'cyclic_shields': cyclic, 'product_group_shields': groups,
        'cyclic_states': states,
        'total_shield_entries_checked': sum(row['shield_entries_checked'] for row in cyclic + groups),
        'negative_control': 'Using the difference shield on both sides fails for Z/3Z.',
    }, indent=2))


if __name__ == '__main__':
    main()
