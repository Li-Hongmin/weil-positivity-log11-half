"""Recompute positivity of the retained small Schur matrices with Arb."""
from pathlib import Path
import numpy as np
from flint import arb, ctx


def check():
    ctx.prec = 256
    root = Path(__file__).resolve().parent / "certificate" / "results"
    paths = sorted(root.glob("*.npz"))
    assert len(paths) == 2
    for path in paths:
        data = np.load(path, allow_pickle=False)
        permutation = [int(i) for i in data["perm"]]
        expected = 7 if str(data["parity"]) == "+" else 6
        assert len(permutation) == expected
        assert sorted(permutation) == list(range(expected))
        raw = data["S"]
        matrix = [[arb(str(raw[i, j])) for j in permutation] for i in permutation]
        for k in range(expected):
            pivot = matrix[k][k]
            assert pivot.lower() > 0, (path, k, pivot)
            for i in range(k + 1, expected):
                for j in range(i, expected):
                    value = matrix[i][j] - matrix[i][k] * matrix[j][k] / pivot
                    matrix[i][j] = value
                    matrix[j][i] = value
        print(path.name, "stored core matrix: PASS; last pivot", pivot)


if __name__ == "__main__":
    check()
