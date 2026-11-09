"""One-shot blockwise contraction of the K=400 safe base.

This avoids the slow path which first parses two dense 400 by 400 string
arrays and only then applies the seven Householder reflectors.  Only the
PP/PQ/QQ blocks are materialised; reflectors act on the first 280 rows.
"""
import time
import sys
import numpy as np
from flint import arb, arb_mat, ctx

BITS = int(sys.argv[sys.argv.index('--precision')+1]) if '--precision' in sys.argv else 192
CACHE_PREFIX = (sys.argv[sys.argv.index('--cache-prefix')+1]
                if '--cache-prefix' in sys.argv else 'log10')
ctx.prec = BITS
R, Q, E = 280, 120, 400
N500_CLEAN = "--n500-clean" in sys.argv
ODD = "--odd" in sys.argv
tag = "odd" if ODD else "even"
parity = "-" if ODD else "+"
core = 6 if ODD else 7
old = R-core
src = np.load(f"scratch/{CACHE_PREFIX}_{tag}_K400_safe_base_arb_{BITS}.npz")
ex_name = ("scratch/log10_even_R280_N500_verified_graph_arb_256.npz"
           if N500_CLEAN and not ODD and CACHE_PREFIX == 'log10' else
           f"scratch/{CACHE_PREFIX}_{tag}_R280_N500_verified_graph_arb_{BITS}.npz"
           if N500_CLEAN and ODD else
           f"scratch/{CACHE_PREFIX}_{tag}_R280_N500_verified_graph_arb_{BITS}.npz"
           if N500_CLEAN else
           "scratch/log10_even_R280_N400_exact_small_map_256.npz")
ex = np.load(ex_name)
whist = ex["whist"] if N500_CLEAN else ex["whist_arb"]
if N500_CLEAN:
    assert (str(src["parity"]) == parity and int(src["R"]) == E
            and int(src["N"]) == 500 and int(src["precision"]) == BITS)
    assert int(ex["N"]) == 500 and int(ex["R"]) == R
t0 = time.time()


def parse_block(raw, i0, i1, j0, j1):
    return arb_mat([[arb(str(raw[i, j])) for j in range(j0, j1)]
                    for i in range(i0, i1)])


def contract(raw, label):
    print(label, "PARSE_PP", flush=True)
    pp = parse_block(raw, 0, R, 0, R)
    print(label, "PARSE_PQ", flush=True)
    pq = parse_block(raw, 0, R, R, E)
    print(label, "PARSE_QQ", flush=True)
    qq = parse_block(raw, R, E, R, E)
    # Congruence by the exact same ordered reflectors used by the certificate.
    for row in whist:
        h = arb_mat([[arb(str(row[i]))] for i in range(R)])
        hh = (h.transpose() * h)[0, 0]
        ph = pp * h
        pp = (pp - h * ((arb(2) / hh) * (h.transpose() * pp))
              - ph * ((arb(2) / hh) * h.transpose())
              + h * ((arb(4) * (h.transpose() * ph)[0, 0] / (hh * hh))
                     * h.transpose()))
        pq = pq - h * ((arb(2) / hh) * (h.transpose() * pq))
    # Fixed-width Unicode avoids the pathological pickle/object-array load
    # time of the historical cache while preserving every Arb ball string.
    out = np.empty((E, E), dtype='<U96')
    for i in range(R):
        for j in range(R): out[i, j] = str(pp[i, j])
        for j in range(Q):
            out[i, R + j] = str(pq[i, j])
            out[R + j, i] = str(pq[i, j])
    for i in range(Q):
        for j in range(Q): out[R + i, R + j] = str(qq[i, j])
    print(label, "DONE_SECONDS", time.time() - t0, flush=True)
    assert max(len(x) for x in out.flat) < 96
    return out


S = contract(src["S"], "S")
Rc = contract(src["Rcov"], "RCOV")
dst = (f"scratch/{CACHE_PREFIX}_{tag}_K400_N500_clean_full_contract_base_arb_{BITS}.npz"
       if N500_CLEAN else
       "scratch/log10_even_K400_full_contract_base_arb_192.npz")
np.savez(dst, S=S, Rcov=Rc,
         source=np.array(f"{CACHE_PREFIX}_{tag}_K400_safe_base_arb_{BITS}.npz"),
         householder_source=np.array(ex_name.split("/")[-1]),
         N=np.array(500 if N500_CLEAN else 400), R=np.array(R),
         clean=np.array(N500_CLEAN), parity=np.array(parity),
         core=np.array(core), old=np.array(old), precision=np.array(BITS),
         cache_prefix=np.array(CACHE_PREFIX), a=np.array(str(src['a'])),
         mu=np.array(str(src['mu'])))
if N500_CLEAN:
    np.savez(f'scratch/{CACHE_PREFIX}_{tag}_K400_N500_clean_scalar_blocks_arb_{BITS}.npz',
             S00=S[:core,:core], S01=S[:core,core:280], S02=S[:core,280:400],
             S11=S[core:280,core:280], S12=S[core:280,280:400],
             S22=S[280:400,280:400],
             source=np.array(dst), N=np.array(500), R=np.array(280),
             K=np.array(400), precision=np.array(BITS), clean=np.array(True),
             parity=np.array(parity), core=np.array(core), old=np.array(old),
             a=np.array(str(src['a'])),
             mu=np.array(str(src['mu'])), cache_prefix=np.array(CACHE_PREFIX))
print("WROTE", dst, "TOTAL_SECONDS", time.time() - t0, flush=True)
