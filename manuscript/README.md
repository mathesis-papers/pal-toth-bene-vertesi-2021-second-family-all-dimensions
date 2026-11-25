# Second-family PPT states in every dimension

This package gives explicit data for the existence assertion in Section V,
paragraph after Eq. (69), p. 023101-8, of K. F. Pál, G. Tóth, E. Bene and T. Vértesi,
*Bound entangled singlet-like states for quantum metrology*, Physical Review
Research **3**, 023101 (2021), DOI: 10.1103/PhysRevResearch.3.023101.

## Contents

- `note.tex`: complete mathematical note with embedded bibliography.
- `verify.py`: standard-library Python 3 exact finite checks.
- `corollary_checks.py`: exact noise-threshold and weighted-Hamiltonian checks.
- `verification.json` and `corollary-verification.json`: finite-check results.
- `metadata.json`: title, authors, license, coverage and fixed scientific source.
- `README.md`: build and calculation instructions.

## Build

Use Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively, use a standard LaTeX installation:

```sh
pdflatex note.tex
pdflatex note.tex
```

The bibliography is embedded, so no separate BibTeX step is needed.

## Exact finite checks

Use Python 3.10 or later, with no third-party dependencies:

```sh
python3 verify.py > verification.json
python3 corollary_checks.py > corollary-verification.json
```

The output reports `status: pass`. Its SHA256 is
`f3a86bba00467831f3166113d7a77c136ad58605ed993a04aceb33d7ad690428`.
It checks 321,592 shield matrix entries, covering cyclic dimensions 2 through
16 and product groups with moduli (2,2), (2,3), (2,2,2), (3,3), and (2,2,2,2).
For every cyclic dimension in that range it also checks the full density
matrix partial-transpose identity, trace, orthogonality, normalization,
Hamiltonian moments and the spectral Fisher-information formula. Four
visibility values per dimension check the white-noise expression. Arithmetic
uses rational coefficients in Q[r]/(r²-d), specialized at r = sqrt(d).
An explicit negative control rejects use of the difference shield on both
sides for the group Z/3Z.

The second program checks the threshold polynomial in a quadratic extension,
90 weighted spectral cases, and 360 weighted white-noise cases for dimensions
2 through 16. It checks six rational pairs of Hamiltonian weights, including
equal, opposite, zero-component and fractional weights. Its scope is finite;
the all-dimension proofs are in the note.

## Mathematical coverage

The all-dimension proof is in `note.tex`. The pinned declaration `result` in
[SecondFamilyAllDimensions.lean](https://github.com/the-omega-institute/trureturing/blob/a2d277ea0980f1230eff2c7e27838372effccb30/D5/S3/Quantum/Entanglement/BoundEntanglement/SecondFamilyAllDimensions.lean)
proves the six conditions defining `IsSecondFamilyData`: orthogonal matrices,
an orthonormal coefficient basis, factorization, multiplicity, orthonormal
shield vectors and the shield partial-transpose identity. The convention for
partial transposition is in
[Bipartite.lean](https://github.com/the-omega-institute/trureturing/blob/a2d277ea0980f1230eff2c7e27838372effccb30/D5/S3/Quantum/Entanglement/AbsolutePPT/Bipartite.lean).
The density-matrix spectra, Fisher information, variance equality, white-noise
formula and abelian-group extension have written proofs in the note.
Corollary 2 gives the exact visibility threshold
`v_* = ((2p-1)+sqrt((2p-1)^2+16p^2))/(8p^2)`, with
`p=sqrt(d)/(1+sqrt(d))`, and limiting noise tolerance `(7-sqrt(17))/8`.
Corollary 3 gives `F_Q(rho,H_c)=4p(c1+c2)^2` and gain at most `2p`
for `H_c=c1 Z_A+c2 Z_B`, with equality exactly at equal weights.
Corollary 4 extends the construction to every finite abelian group,
gives product decompositions of the normalized shield projectors, and proves
that same-shield partial-transpose invariance holds exactly for exponent-two groups.
The Python checks cover only the displayed finite ranges.
Conjecture 1 of Pál et al., concerning global optimality among all PPT states,
is not covered by the construction.

## Authors and license

Haobo Ma — ChronoAI Pte Ltd, Singapore.

Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore.

Copyright 2026 Haobo Ma and Wenlin Zhang.
The note and accompanying original materials are licensed under
[Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/),
identifier `cc-by-4.0`.
