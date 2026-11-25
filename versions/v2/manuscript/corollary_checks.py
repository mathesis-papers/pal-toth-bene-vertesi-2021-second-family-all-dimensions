#!/usr/bin/env python3
"""Exact finite checks of the noise threshold and weighted Hamiltonian formulas."""

from fractions import Fraction
import json

from verify import Quadratic


def inverse(value):
    norm = value.rational ** 2 - value.dimension * value.radical ** 2
    assert norm != 0
    return Quadratic(value.dimension, value.rational / norm, -value.radical / norm)


def pair_add(left, right):
    return left[0] + right[0], left[1] + right[1]


def pair_multiply(left, right, discriminant):
    return (left[0] * right[0] + discriminant * left[1] * right[1],
            left[0] * right[1] + left[1] * right[0])


def check_dimension(dimension):
    weight_p = Quadratic(dimension, Fraction(dimension, dimension - 1),
                         Fraction(-1, dimension - 1))
    weight_q = 1 - weight_p
    zero = Quadratic(dimension)
    half = Quadratic(dimension, Fraction(1, 2))
    linear = weight_p - half
    discriminant = linear * linear + 4 * weight_p * weight_p
    denominator_inverse = inverse(4 * weight_p * weight_p)
    threshold = linear * denominator_inverse, denominator_inverse
    square = pair_multiply(threshold, threshold, discriminant)
    polynomial = pair_add(pair_multiply((2 * weight_p * weight_p, zero),
                                       square, discriminant),
                          pair_multiply((-linear, zero), threshold, discriminant))
    polynomial = pair_add(polynomial, (-half, zero))
    assert polynomial == (zero, zero)
    weight_pairs = [(1, 1), (1, 2), (1, -1), (0, 1), (2, -3),
                    (Fraction(2, 3), Fraction(5, 7))]
    visibilities = [Fraction(0), Fraction(1, 5), Fraction(1, 2), Fraction(1)]
    weighted_cases = 0
    noise_cases = 0
    for first, second in weight_pairs:
        eigenvalue = weight_p * Fraction(1, dimension ** 2)
        pair_element = Fraction((first + second) ** 2)
        spectral_sum = zero
        second_moment = zero
        mean = zero
        for row in range(dimension):
            for column in range(dimension):
                plus_support = [(0, 0, row, column),
                                (1, 1, column, (row + column) % dimension)]
                energies = [first * (1 - 2 * basis[0])
                            + second * (1 - 2 * basis[1]) for basis in plus_support]
                diagonal = Fraction(sum(energies), 2)
                off_diagonal = Fraction(energies[0] - energies[1], 2)
                assert diagonal == 0
                assert off_diagonal ** 2 == pair_element
                spectral_sum += 4 * eigenvalue * off_diagonal ** 2
                mean += eigenvalue * diagonal
                second_moment += eigenvalue * Fraction(sum(energy ** 2 for energy in energies), 2)
        for key_a, key_b in [(0, 1), (1, 0)]:
            energy = first * (1 - 2 * key_a) + second * (1 - 2 * key_b)
            mean += weight_q * Fraction(energy, 2)
            second_moment += weight_q * Fraction(energy ** 2, 2)
        assert mean == 0
        assert spectral_sum == 4 * weight_p * (first + second) ** 2
        assert 4 * second_moment - spectral_sum == 4 * weight_q * (first - second) ** 2
        assert (first + second) ** 2 + (first - second) ** 2 == 2 * (first ** 2 + second ** 2)
        weighted_cases += 1
        for visibility in visibilities:
            floor = Fraction(1 - visibility, 4 * dimension ** 2)
            plus_eigenvalue = visibility * eigenvalue + floor
            minus_eigenvalue = Quadratic(dimension, floor)
            difference = plus_eigenvalue - minus_eigenvalue
            spectral_numerator = 4 * dimension ** 2 * difference * difference * pair_element
            spectral_denominator = plus_eigenvalue + minus_eigenvalue
            target_numerator = spectral_sum * 2 * weight_p * visibility ** 2
            target_denominator = 1 + (2 * weight_p - 1) * visibility
            assert spectral_numerator * target_denominator == target_numerator * spectral_denominator
            noise_cases += 1
    return {'dimension': dimension, 'threshold_polynomial_zero': True,
            'weighted_spectral_cases': weighted_cases, 'weighted_noise_cases': noise_cases,
            'variance_defect_identity': True, 'equal_weight_gain_maximum': True}


def main():
    dimensions = [check_dimension(dimension) for dimension in range(2, 17)]
    print(json.dumps({
        'status': 'pass',
        'arithmetic': 'Exact rational coefficients in Q[r]/(r^2-d), extended by s^2=(p-1/2)^2+4p^2 for the threshold.',
        'coverage': 'Cyclic dimensions 2 through 16; six weight pairs and four visibilities per pair. Universal claims have written proofs in note.tex.',
        'dimensions': dimensions,
        'total_weighted_spectral_cases': sum(row['weighted_spectral_cases'] for row in dimensions),
        'total_weighted_noise_cases': sum(row['weighted_noise_cases'] for row in dimensions),
        'threshold_polynomial_cases': len(dimensions),
    }, indent=2))


if __name__ == '__main__':
    main()
