# log11 endpoint certificate, both parity sectors (run of 10 September 2026)

Fixed parameters: a = log(11)/2, K = 400, N = 500, R = 280, response family
F = [500, 850), q_f = 8000, Arb precision 256 bits.

For each parity the run starts from an empty `scratch/` directory, generates the
safe base, the verified Householder graph, the clean contracted base and the
complete response-family Parseval data, and then runs the sharp checker. The large
parent files stay in the unversioned `scratch/` directory; their SHA-256 values are
listed in `PARENT-SHA256.txt`. This bundle contains the construction sources, the
stage logs and the final small interval results.

The checker forms
`G = P0 + RP/t400 + BF*BF.T - G1 - (Rtot - R1)/t8000`, where
`Rtot = RP + XF*BF.T + BF*XF.T + BF*GQ*BF.T`, and builds `U_q = c_q + BF*v_q`
term by term. The generator keeps the negative odd pole term; the family exporter
stores all cross terms and all `gram_q`. The trigamma series shows that `t_q` is
strictly increasing, and both computed values of `t400` are strictly positive.
Since `Rtot - R1` is the Parseval Gram matrix of the remaining `U_q U_q^*`, the
tail weight has the correct direction.

Even sector: the minimum pivot of the 393-dimensional complementary block has
lower endpoint about 0.04079956, and the last core pivot has lower endpoint about
3.6824e-50. Odd sector: the corresponding values are about 0.00733669 for the
394-dimensional block and 3.6435e-45 for the last core pivot. All pivots are
strictly positive. The bundle certifies the endpoint (log 11)/2 in both parity
sectors; it does not prove RH.

The supplied checker is a later revision of the one that produced the stored
logs. It requires exactly seven even or six odd core pivots, each with a strictly
positive lower endpoint, and exits with nonzero status on a negative or
zero-containing pivot. For this reason the stored checker logs contain one
diagnostic line, `FROZEN_NEGATIVE_WITNESS_SAVED None`, which the supplied checker
prints only when a negative witness is requested. `verify_frozen_results.py`
checks the parameters, the pivot counts, the strict signs and the bundle hashes;
it compares the manifest with the fixed set of 21 expected relative paths before
hashing.
