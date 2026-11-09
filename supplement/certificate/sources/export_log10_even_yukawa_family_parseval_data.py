"""Export exact Parseval data for a finite family of full even Yukawa atoms.

This file is executed inside an endpoint contract-cache generator
after the clean scalar K400/N500 objects ``Bpp``, ``HS`` and ``mu`` exist.
The driver must set ``YUKAWA_EXPORT_LO`` and ``YUKAWA_EXPORT_HI``.

No complement truncation is used.  If U is the full atom column family,
then the exported matrices satisfy

    b       = P U,
    h       = P B0 U,
    x_v     = P B0 Q (Q U) = h - (P B0 P)b,
    gram_q  = (Q U)^*(Q U) = U^*U - b^*b,
    e       = b - x_v/mu,
    metric  = I + gram_q/mu.

Consequently the exact scalar-complement Woodbury gain is
``e * metric^{-1} * e^*``.  All cross terms between distinct atoms are in
``gram_full``/``metric``; treating atoms independently would be wrong.
"""
import numpy as np

assert PAR in {'+', '-'} and R == 400 and N == 500
tag = 'even' if PAR == '+' else 'odd'
cache_prefix = str(globals().get('CACHE_PREFIX', 'log10'))
lo = int(YUKAWA_EXPORT_LO)
hi = int(YUKAWA_EXPORT_HI)
assert 500 <= lo < hi
count = hi - lo

atoms = [rankrep(n) for n in range(lo, hi)]

# P U.
b = arb_mat([[inner(BAS[i], atoms[j], arb(0), a)
              for j in range(count)] for i in range(R)])

# P times the prime/reflection part of B0 U.  TB already contains the
# transformed physical P-basis rows of every retained prime term.
hp = arb_mat(R, count)
for coeff, xlo, xhi, rb in TB:
    for i in range(R):
        for j in range(count):
            hp[i, j] += coeff * inner(rb[i], atoms[j], xlo, xhi)

# Inner products with all retained n<500 Yukawa columns and the polar column.
q = arb_mat([[inner(FUN[k], atoms[j], arb(0), a)
              for j in range(count)] for k in range(M)])
h = hp + HS * q
x_v = h - Bpp * b
e = b - x_v / mu

# Exact full-space and Parseval-complement atom Grams.
gram_full = arb_mat([[inner(atoms[i], atoms[j], arb(0), a)
                     for j in range(count)] for i in range(count)])
gram_q = gram_full - b.transpose() * b
metric = arb_mat(count, count)
for i in range(count):
    for j in range(count):
        metric[i, j] = gram_q[i, j] / mu + (arb(1) if i == j else arb(0))
metric = (metric + metric.transpose()) / 2

def string_matrix(mat):
    # Fixed-width unicode avoids the very slow object-array deserialization
    # that affected the historical full-contract cache.
    vals = [[str(mat[i, j]) for j in range(mat.ncols())]
            for i in range(mat.nrows())]
    width = max((len(x) for row in vals for x in row), default=1)
    return np.asarray(vals, dtype=f'<U{width}')

out = (f'scratch/{cache_prefix}_{tag}_K400_N500_yukawa_{lo}_{hi}_'
       f'parseval_arb_{ctx.prec}.npz')
np.savez(
    out,
    b=string_matrix(b),
    h=string_matrix(h),
    cross=string_matrix(x_v),
    e=string_matrix(e),
    gram_full=string_matrix(gram_full),
    gram_q=string_matrix(gram_q),
    metric=string_matrix(metric),
    t400=np.array(str(mu)),
    a=np.array(str(a)),
    source_model=np.array(f'{cache_prefix}_{tag}_K400_N500_scalar_lift'),
    cache_prefix=np.array(cache_prefix),
    parity=np.array(PAR), R=np.array(R), N=np.array(N),
    atom_lo=np.array(lo), atom_hi=np.array(hi),
    precision=np.array(ctx.prec),
    exact_parseval=np.array(True),
)
print('EXPORTED_YUKAWA_FAMILY_PARSEVAL', out, flush=True)
print('YUKAWA_ATOM_RANGE', lo, hi, 'COUNT', count, flush=True)
raise SystemExit(0)
