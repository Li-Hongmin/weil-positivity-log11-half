"""Grid-free four-step Schur certificate for log(11)/2 or log(13)/2.

Endpoints are represented as exact
Q-linear combinations of logarithms of distinct primes.  Hence equality of
the tuples is exact (unique factorisation), while Arb certifies their order.
"""
from fractions import Fraction as F
import sys

from flint import arb, ctx


ctx.prec = int(sys.argv[1]) if len(sys.argv) > 1 else 256
cross_endpoint = "--n11-cross" in sys.argv
max_degree = int(sys.argv[sys.argv.index("--max-degree") + 1]) if "--max-degree" in sys.argv else 4
if max_degree < 1:
    raise SystemExit("--max-degree must be positive")
endpoint_prime = int(sys.argv[2]) if len(sys.argv) > 2 else 11
if cross_endpoint:
    # b_cross = log(11) - (1/2)log(201/20)
    # = log(11) + log(2) + (1/2)log(5) - (1/2)log(3) - (1/2)log(67).
    endpoint_prime = 13  # activate the n=11 delay list below
if endpoint_prime not in (11, 13):
    raise SystemExit("endpoint prime must be 11 or 13")

PRIMES = (2, 3, 5, 7, 11, 13, 67) if cross_endpoint else (2, 3, 5, 7, 11, 13)
ZERO = (F(0),) * len(PRIMES)
NAMES = ("x", "y", "z", "w", "v", "r", "s", "q")


def vec(*entries):
    return tuple(F(x) for x in entries)


def add(left, right):
    return tuple(x + y for x, y in zip(left, right))


def neg(value):
    return tuple(-x for x in value)


def sub(left, right):
    return add(left, neg(right))


def mul(scalar, value):
    return tuple(F(scalar) * x for x in value)


def evaluate(value):
    return sum(
        (
            arb(coefficient.numerator)
            / coefficient.denominator
            * arb(prime).log()
            for coefficient, prime in zip(value, PRIMES)
        ),
        arb(0),
    )


def less(left, right):
    # exact equality is decided on the rational vectors; every other
    # comparison must be decided by Arb, otherwise abort.
    if left == right:
        return False
    lv, rv = evaluate(left), evaluate(right)
    if lv < rv:
        return True
    if lv > rv:
        return False
    raise ArithmeticError(("undecided Arb comparison", left, right, lv, rv))


def inside(value, lower, upper):
    return less(lower, value) and less(value, upper)


L2 = vec(1, 0, 0, 0, 0, 0)
L3 = vec(0, 1, 0, 0, 0, 0)
L5 = vec(0, 0, 1, 0, 0, 0)
L7 = vec(0, 0, 0, 1, 0, 0)
L11 = vec(0, 0, 0, 0, 1, 0)
L13 = vec(0, 0, 0, 0, 0, 1)
A = (
    vec(1, F(-1, 2), F(1, 2), 0, 1, 0, F(-1, 2))
    if cross_endpoint
    else mul(F(1, 2), L11 if endpoint_prime == 11 else L13)
)

DELAYS = {
    2: L2,
    3: L3,
    4: mul(2, L2),
    5: L5,
    7: L7,
    8: mul(3, L2),
    9: mul(2, L3),
}
if endpoint_prime == 13:
    DELAYS[11] = L11
WEIGHT_NAME = {2: "x", 3: "y", 4: "z", 5: "w", 7: "v", 8: "r", 9: "s"}
if endpoint_prime == 13:
    WEIGHT_NAME[11] = "q"

# An edge is (weight index, affine sign, affine constant, support lo, support hi, label).
EDGES = []
for n, delay in DELAYS.items():
    if less(delay, A):
        EDGES.extend(
            [
                (NAMES.index(WEIGHT_NAME[n]), 1, delay, ZERO, sub(A, delay), f"S{n}+"),
                (NAMES.index(WEIGHT_NAME[n]), 1, neg(delay), delay, A, f"S{n}-"),
            ]
        )
for n, delay in DELAYS.items():
    lower = ZERO if less(delay, A) else sub(delay, A)
    upper = delay if less(delay, A) else A
    EDGES.append(
        (NAMES.index(WEIGHT_NAME[n]), -1, delay, lower, upper, f"J{n}")
    )


def image(edge, value):
    return add(mul(edge[1], value), edge[2])


def inverse_image(edge, value):
    return mul(edge[1], sub(value, edge[2]))


def poly_add(target, source):
    for monomial, coefficient in source.items():
        target[monomial] = target.get(monomial, 0) + coefficient


def poly_times_variable(source, variable):
    result = {}
    for monomial, coefficient in source.items():
        exponent = list(monomial)
        exponent[variable] += 1
        exponent = tuple(exponent)
        result[exponent] = result.get(exponent, 0) + coefficient
    return result


def walk_polynomial(start, steps):
    zero_monomial = (0,) * len(NAMES)
    states = {start: {zero_monomial: 1}}
    for _ in range(steps):
        following = {}
        for value, polynomial in states.items():
            for edge in EDGES:
                if inside(value, edge[3], edge[4]):
                    destination = image(edge, value)
                    poly_add(
                        following.setdefault(destination, {}),
                        poly_times_variable(polynomial, edge[0]),
                    )
        states = following
    total = {}
    for polynomial in states.values():
        poly_add(total, polynomial)
    return total


WEIGHTS = (
    arb(2).log() / arb(2).sqrt(),
    arb(3).log() / arb(3).sqrt(),
    arb(2).log() / 2,
    arb(5).log() / arb(5).sqrt(),
    arb(7).log() / arb(7).sqrt(),
    arb(2).log() / arb(8).sqrt(),
    arb(3).log() / 3,
    arb(11).log() / arb(11).sqrt(),
)


def polynomial_value(polynomial):
    total = arb(0)
    for exponent, coefficient in polynomial.items():
        term = arb(coefficient)
        for weight, power in zip(WEIGHTS, exponent):
            term *= weight**power
        total += term
    return total


# D_4 can change only at a support endpoint or at one of its inverse images
# under a word of length at most three.
cuts = {ZERO, A}
for edge in EDGES:
    cuts.update((edge[3], edge[4]))
for _ in range(max_degree - 1):
    old_cuts = tuple(cuts)
    for edge in EDGES:
        for breakpoint in old_cuts:
            candidate = inverse_image(edge, breakpoint)
            if candidate in (ZERO, A) or inside(candidate, edge[3], edge[4]):
                cuts.add(candidate)

cuts = sorted(cuts, key=lambda value: float(evaluate(value).mid()))
for lower, upper in zip(cuts, cuts[1:]):
    assert less(lower, upper), (evaluate(lower), evaluate(upper))

old_rho = arb("2.978071" if endpoint_prime == 11 else "2.931158")
print("bits", ctx.prec)
print("endpoint_prime", endpoint_prime)
print("max_degree", max_degree)
print("active_n", sorted(DELAYS), "edge_count", len(EDGES))
print("cuts", len(cuts), "cells", len(cuts) - 1)
for steps in range(1, max_degree + 1):
    rows = []
    for lower, upper in zip(cuts, cuts[1:]):
        midpoint = mul(F(1, 2), add(lower, upper))
        polynomial = walk_polynomial(midpoint, steps)
        rows.append((lower, upper, polynomial_value(polynomial), polynomial))

    best = max(rows, key=lambda row: float(row[2].mid()))
    supremum = best[2]
    for lower, upper, value, polynomial in rows:
        assert polynomial == best[3] or value < supremum, (
            evaluate(lower),
            evaluate(upper),
            value,
            supremum,
        )

    root = supremum ** (arb(1) / steps)
    print(f"sup_D{steps}", supremum)
    print(f"root_D{steps}", root, "below_old_rho", root < old_rho)
    if steps == 4:
        print("max_cell_D4", evaluate(best[0]), evaluate(best[1]))
        print("max_polynomial_terms_D4", len(best[3]))
        target = arb("2.931158")
        assert root < target, (root, target)
        print("root_D4_below_2.931158", root < target, "margin_lower", (target - root).lower())
