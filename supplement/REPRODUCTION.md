# Reproducing the endpoint certificate

The manuscript proves the comparison with the complete Weil form. The files
in `certificate/` are the existing fixed endpoint certificate and original
construction sources. Their historical names contain `log5`, `log8`, or
`log10`; the commands below fix the endpoint to 11 and the precision to 256.
No private repository, unpublished parent matrix, or external data is needed.

## Short checks

From the extracted supplement root:

```sh
uv run --with numpy python certificate/verify_frozen_results.py
uv run --with numpy --with python-flint python verify_retained_core_matrices.py
uv run --with python-flint python certify_log11_prime_four_step_degree.py 256 11
uv run --with sympy python verify_green_response_identities.py
```

Expected results: both frozen parity results pass metadata, complete pivot
counts, signs and their existing file bindings; the graph calculation has
200 cells and fourth-root bound below 2.931158; both Green identities and
the digamma summand normalization pass exact symbolic checks.
The additional core check recomputes interval elimination on the actual
stored seven- and six-dimensional Schur matrices, rather than reading only
their retained pivot lists.

The frozen verifier checks retained results, not the full generation of
their matrices. The graph enumerates exact logarithmic endpoints and uses
outward Arb comparisons. The symbolic check validates algebraic identities,
not form-domain or infinite-dimensional assertions.

## Full sequential regeneration

Create a new directory so that the frozen certificate remains untouched:

```sh
mkdir -p work/computations work/scratch
cp certificate/sources/*.py work/computations/
cp certificate/sources/export_log10_even_yukawa_family_parseval_data.py work/scratch/
cd work
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
```

Run each command below sequentially. The first command creates the full
safe base. The second creates verified Householder coordinates. The third
contracts the base in those coordinates. The fourth evaluates the entire
retained physical response family, including all cross Grams. The fifth
checks the complete Fourier-complement lower model using the weighted tail.

Even sector:

```sh
uv run --with numpy --with python-flint python computations/generate_log10_endpoint_contract_cache.py 256 400 500 + --endpoint 11 --export-full-safe-matrix
uv run --with numpy --with python-flint python computations/generate_log10_endpoint_contract_cache.py 256 280 500 + --endpoint 11 --export-verified-graph
uv run --with numpy --with python-flint python computations/extract_log10_even_K400_full_contract_base.py --cache-prefix log11 --precision 256 --n500-clean
uv run --with numpy --with python-flint python computations/generate_log10_endpoint_contract_cache.py 256 400 500 + --endpoint 11 --export-yukawa-family 500 850
uv run --with numpy --with python-flint python computations/log10_even_K400_sharp_safe_arb.py --full400 --clean-no-yukawa --family scratch/log11_even_K400_N500_yukawa_500_850_parseval_arb_256.npz 350 --qf 8000 --cache-prefix log11 256 11 1
```

Odd sector:

```sh
uv run --with numpy --with python-flint python computations/generate_log10_endpoint_contract_cache.py 256 400 500 - --endpoint 11 --export-full-safe-matrix
uv run --with numpy --with python-flint python computations/generate_log10_endpoint_contract_cache.py 256 280 500 - --endpoint 11 --export-verified-graph
uv run --with numpy --with python-flint python computations/extract_log10_even_K400_full_contract_base.py --cache-prefix log11 --precision 256 --n500-clean --odd
uv run --with numpy --with python-flint python computations/generate_log10_endpoint_contract_cache.py 256 400 500 - --endpoint 11 --export-yukawa-family 500 850
uv run --with numpy --with python-flint python computations/log10_even_K400_sharp_safe_arb.py --full400 --clean-no-yukawa --odd --family scratch/log11_odd_K400_N500_yukawa_500_850_parseval_arb_256.npz 350 --qf 8000 --cache-prefix log11 256 11 1
```

Each final checker must exit successfully and report all seven even or six
odd core pivots strictly positive. The larger complementary blocks must
also pass their complete interval Cholesky checks. Individual midpoint
eigenvalues are not certificates.

The old runs used Python 3.10.12 on a dedicated compute host. Exact historical
NumPy and python-flint version numbers were not retained in this small
bundle. The submission includes the versions used for its short checks in
`check-environment.txt`; those short checks do not establish a complete
regeneration with those versions. Different floating-point libraries can
select different Householder proposals. Their validity must still be checked
by the interval stages, rather than requiring bit-for-bit identical proposals.

## Resources and evidence scope

Use one process and one parity at a time. The retained even sharp-prefix
stage took approximately 1474 seconds before the final complement solve;
base construction and the other stages add substantial time. The large
string-valued interval parent matrices are generated under `work/scratch/`
and are not bundled. A full regeneration was not run during this submission
preparation, so no new runtime or peak-memory measurement is claimed.

`prime-graph-check.txt` and `short-checks.txt` record the short checks actually
run for this submission. The files within `certificate/results/` are the
historical full-run records. The theorem's infinite-dimensional assertions
depend on the manuscript's analytic proof as well as the finite certificate.
