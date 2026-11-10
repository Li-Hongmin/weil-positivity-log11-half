# Cutoff-free positivity of localized Weil forms through the (log 11)/2 threshold

[![Verify certificates](https://github.com/Li-Hongmin/weil-positivity-log11-half/actions/workflows/verify.yml/badge.svg)](https://github.com/Li-Hongmin/weil-positivity-log11-half/actions/workflows/verify.yml)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.24054983-blue)](https://doi.org/10.5281/zenodo.24054983)

**An AI-generated manuscript with a computer-assisted certificate, released for human verification.**

Hongmin Li (Institute of Science Tokyo; The University of Tokyo) · lihongmin@edu.k.u-tokyo.ac.jp

This repository holds the manuscript in the version accepted on rStation
(version 2, 9 November 2026: <https://rstation.ai/problems/riemann_hypothesis/papers/45/v2>;
archived at Zenodo, DOI [10.5281/zenodo.24076153](https://doi.org/10.5281/zenodo.24076153)), the supplement
with its verification programs and stored certificate data, and a
continuous-integration workflow that reruns the short checks on a clean
machine.

## What is claimed

Weil's criterion expresses the Riemann hypothesis as positivity of an explicit
quadratic form on every compact support. The manuscript studies the complete
geometric Weil form on complex-valued functions supported in $(-a,a)$, in the
normalization of its equation (2.1), and its Friedrichs operator $A_a$.

> **Theorem 1.1.** Let $A_a$ be the Friedrichs operator of the geometric Weil
> form on $(-a,a)$ in the normalization (2.1). Then
> $\inf\sigma(A_a)>0$ for $0<a\le a_*=\tfrac12\log 11=1.198947\ldots$.
> Both parity restrictions are strictly positive. If the first degeneracy
> radius is defined by $a_{\mathrm{deg}}=\inf\{a>0:\inf\sigma(A_a)\le0\}$,
> with the infimum of the empty set taken to be infinity, then
> $a_{\mathrm{deg}}>a_*$.

Method, in one sentence: the archimedean part is written as a diagonal
digamma operator plus positive rank-one boundary responses of exponential
Green kernels; a finite response family is retained, its cross terms with the
entire Fourier complement are evaluated by Parseval and bounded by a monotone
diagonal weight, and a 400-dimensional Schur matrix is then certified positive
with 256-bit interval arithmetic in each parity sector.

The manuscript gives a proof of Theorem 1.1. That proof has not yet been
checked by a human expert.

**Not claimed:**

- not the Riemann hypothesis, a new zero-free region, or all-scale Weil
  positivity;
- nothing for any support radius $a>\tfrac12\log 11$; a larger finite
  response family has no automatic implication for larger radii;
- no explicit coercivity constant and no spectral-gap value: the small
  elimination pivots reported in the manuscript are not eigenvalue bounds;
- no judgement on claims by other authors that have not been verified here
  (see "Relation to prior work").

## How this manuscript was produced

The manuscript, the certificate construction and the verification programs
were produced by AI systems under my direction. My role was to set goals,
allocate compute, and decide when to stop. **I am not a specialist in analytic
number theory or operator theory, and I cannot independently vouch for the
correctness of the proof.**

I am releasing it for exactly that reason. If you work on Weil positivity,
explicit formulae, or certified numerics, I would be grateful for your
reading, whether the outcome is "correct", "fixable", or "wrong at step X".
I will publicly record whichever it turns out to be.

Substantive verification or correction will be acknowledged. Co-authorship on
any revised version is open to anyone who materially contributes.

## Relation to prior work

The manuscript positions itself as follows; the wording matches its
introduction.

- **Suzuki, arXiv:2606.09096.** Provides the operator realization of the
  form, continuity of the lowest eigenvalue in the support radius, and
  positivity for sufficiently small $a$. The manuscript uses these results
  as inputs and proves nothing about them.
- **Yoshida, Adv. Stud. Pure Math. 21 (1992).** Proved strict positivity of
  the localized form, pole term included, at support radius $\tfrac12\log 2$
  by a finite-matrix-plus-tail argument with computer-assisted interval
  bounds. The finite-matrix-plus-tail strategy goes back to this work.
- **Zhu, arXiv:2608.24827.** Proves an explicit coercivity bound for all
  complex test functions supported in $[-0.8,0.8]$.
- **Liu, certified-weil-positivity (2026).** Gives computer-assisted
  coercivity constants at half-widths $1$ and $17/16$, using source-exact
  block Schur and tail compensation.
- **Desogus, arXiv:2609.20367v3.** Proposes an all-scale argument in the
  real odd channel. Its computer-assisted input is a certificate for the odd
  localized operator at support radius $\tfrac12\log 7$; an arithmetic Schur
  induction over the integer endpoints is then asserted to give odd
  positivity at every support radius and hence RH. That induction has not
  been verified here, and nothing in this manuscript depends on it.

The endpoint $\tfrac12\log 11=1.1989\ldots$ of the present manuscript exceeds
the half-widths $0.8$, $1$ and $17/16$ above. Unlike Zhu and Liu, Theorem 1.1
states no explicit coercivity constant. Block Schur complementation, Parseval's
identity and monotone tail bounds are standard; the manuscript makes no
priority claim for them. Expert opinion on how the present certificate
relates to these works is welcome.

## What the computer verification does and does not cover

The continuous-integration workflow (`.github/workflows/verify.yml`) runs the
four short checks below on every push, on a clean Ubuntu runner with pinned
dependencies (`requirements.txt`). Together they take well under a minute.
Before them it verifies `supplement/CONTENTS.sha256`, the checksum manifest
of the supplement exactly as submitted.

| Script (run from `supplement/`) | What it verifies | Manuscript reference |
|---|---|---|
| `certificate/verify_frozen_results.py` | Loads the two stored endpoint results (even and odd). Checks parameters ($a=\tfrac12\log 11$, $K=400$, $N=500$, family $[500,850)$, $q_f=8000$, 256 bits), the complete pivot lists (seven even, six odd), the stored permutations, strictly positive lower endpoints of every pivot and of the complement minimum, $t_{400}>0$, the monotone tail endpoints, and the SHA-256 of all 21 files in `certificate/`. | Section 7, display (7.1) |
| `verify_retained_core_matrices.py` | Re-runs the symmetric pivoted interval elimination on the stored $7\times7$ and $6\times6$ Schur matrices at 256 bits, rather than only reading the stored pivots. | Section 7 |
| `verify_green_response_identities.py` | Exact symbolic (SymPy) check of the Green boundary-response identity in both parities and of the rational identity behind the digamma summand. | Section 3, displays (3.3) and (3.4) |
| `certify_log11_prime_four_step_degree.py 256 11` | Exhaustive four-step enumeration of the prime operator on 200 cells with exact rational endpoints and Arb ordering; asserts the fourth root is below $2.931158$. CI compares the first 15 output lines with the stored `prime-graph-check.txt`. | Section 6, display (6.2) and bound (4.3) |

**Not reproduced in CI, and not formalized.** This is where a human reader is
needed.

1. The construction of the $400\times400$ matrix $G$ of display (5.4) in each
   parity, including the full safe base, the Householder congruence, the
   Parseval covariance $R_{\mathrm{tot}}$ and the interval Cholesky of the
   393- and 394-dimensional blocks. The repository contains the construction
   sources, the stage logs of the retained run of 10 September 2026, and the
   logs of a complete regeneration on 8 November 2026
   (`supplement/rerun-20261108/`). The regeneration takes roughly an hour
   per parity and produces parent matrices of about 0.7 GB; it is not run in
   CI. Instructions are in `supplement/REPRODUCTION.md`.
2. The identification of the computed matrix with display (5.4), that is,
   that the construction sources implement the operators defined in the
   manuscript. The short checks verify stored outputs, not that the sources
   compute the right objects.
3. Every analytic step: Lemma 3.1 and its convergence argument, the lower
   operator and the passage to the complete form domain in Section 4
   (displays (4.7) to (4.9)), the Parseval and monotone-tail inequalities
   (5.5), the Schur-test reduction (6.1), the zero-extension and min-max
   argument for $a<a_*$, and the use of continuity for $a_{\mathrm{deg}}>a_*$.
4. The cited inputs from Suzuki and from Connes, Consani and Moscovici
   (operator realization, form core, discreteness of the spectrum,
   continuity of the lowest eigenvalue).

The checks rule out arithmetic slips in the retained finite data and confirm
two algebraic identities. They do not constitute a proof of Theorem 1.1.

## Reproducing

Short checks (about 15 seconds in total on a laptop):

```sh
cd supplement
python -m pip install -r ../requirements.txt
shasum -a 256 -c CONTENTS.sha256
python certificate/verify_frozen_results.py
python verify_retained_core_matrices.py
python verify_green_response_identities.py
python certify_log11_prime_four_step_degree.py 256 11 | diff - prime-graph-check.txt
```

Or fork this repository: the GitHub Actions workflow runs the same commands
on a clean machine.

Full regeneration of the endpoint certificate from an empty directory follows
`supplement/REPRODUCTION.md` (NumPy and python-flint only; one parity at a
time; about 0.7 GB of intermediate interval matrices). Different
floating-point environments select different Householder proposals, so the
regenerated pivots differ from the stored ones and must pass the interval
stages on their own; the logs in `supplement/rerun-20261108/` show one such
regeneration.

The PDF is typeset from `paper/paper.md` with `paper/build-pdf.sh`
(pandoc and tectonic); the title and abstract come from `paper/metadata.json`.

## Contents

| Path | What it is |
|---|---|
| `paper/paper.md` | The manuscript as submitted (rStation version 2) |
| `paper/paper.pdf` | The same manuscript typeset (9 pages) |
| `paper/metadata.json` | Submission title, abstract, tags and change summary |
| `supplement/CONTENTS.sha256` | SHA-256 manifest of every supplement file as submitted |
| `supplement/REPRODUCTION.md` | Short checks and full sequential regeneration |
| `supplement/certificate/` | Stored endpoint results, stage logs, construction sources, checksums, frozen-result verifier |
| `supplement/rerun-20261108/` | Logs, timings and checksums of the complete regeneration of 8 November 2026 |
| `supplement/*.py` | Prime enumeration, Green-identity check, retained-core re-elimination |
| `supplement/prime-graph-check.txt`, `short-checks.txt`, `check-environment.txt` | Recorded outputs of the short checks at submission |
| `.github/workflows/verify.yml` | CI: reruns the four short checks with pinned dependencies |

## Cite

Hongmin Li, *Cutoff-free positivity of localized Weil forms through the
(log 11)/2 threshold*, rStation, riemann_hypothesis paper 45, 2026.
<https://rstation.ai/problems/riemann_hypothesis/papers/45>.
Archived at Zenodo, concept DOI
[10.5281/zenodo.24054983](https://doi.org/10.5281/zenodo.24054983).
Source, certificates and verification programs: this repository.

## License

Manuscript and prose: CC BY 4.0. Code and certificate data: MIT. See `LICENSE`.
