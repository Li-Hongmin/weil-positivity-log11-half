"""Verify fixed log11 round-1 artifacts without relying on ``assert``."""
from pathlib import Path
from decimal import Decimal
from fractions import Fraction
import hashlib
import re
import sys
import numpy as np

ROOT = Path(__file__).resolve().parent


def fail(message):
    print(f"FAIL: {message}", file=sys.stderr)
    raise SystemExit(1)


def require(condition, message):
    if not condition:
        fail(message)


def lower_endpoint(value):
    text = str(value).strip()
    match = re.fullmatch(r"\[?([^\[\]]+?) \+/- ([^\[\]]+)\]?", text)
    if match:
        return Fraction(match.group(1)) - Fraction(match.group(2))
    return Fraction(text)


expected = {
    "+": (7, 393, [6, 2, 3, 4, 5, 1, 0]),
    "-": (6, 394, [3, 5, 4, 1, 0, 2]),
}
for parity, (core, tail, expected_perm) in expected.items():
    tag = "even" if parity == "+" else "odd"
    path = ROOT / "results" / (
        f"log11_{tag}_K400_N500_clean_sharp_family_500_850_"
        f"full{tail}_q8000_arb_256.npz"
    )
    try:
        data = np.load(path, allow_pickle=False)
        require(str(data["parity"]) == parity, f"{tag}: parity")
        require(str(data["cache_prefix"]) == "log11", f"{tag}: cache prefix")
        require(Decimal(str(data["endpoint_value"])) == Decimal(11), f"{tag}: endpoint")
        require(bool(data["clean"]), f"{tag}: clean marker")
        require(int(data["N"]) == 500 and int(data["precision"]) == 256,
                f"{tag}: N/precision")
        require(int(data["qf"]) == 8000 and int(data["family_m"]) == 350,
                f"{tag}: qf/family")
        pivots = data["pivots"]
        perm = [int(x) for x in data["perm"]]
        require(len(pivots) == core, f"{tag}: pivot count")
        require(perm == expected_perm, f"{tag}: frozen permutation")
        require(sorted(perm) == list(range(core)), f"{tag}: permutation validity")
        require(lower_endpoint(data["Hmin"]) > 0, f"{tag}: complement lower bound")
        require(lower_endpoint(data["t400"]) > 0, f"{tag}: tail t400 lower bound")
        require(lower_endpoint(data["tqf"]) >= lower_endpoint(data["t400"]),
                f"{tag}: monotone tail endpoints")
        require(all(lower_endpoint(pivot) > 0 for pivot in pivots),
                f"{tag}: strict core pivots")
    except (KeyError, ValueError, TypeError, OSError) as exc:
        fail(f"{tag}: malformed result: {exc}")

    graph_log = (ROOT / "results" / f"{tag}-graph.log").read_text()
    require(graph_log.splitlines()[0].strip() == "tail_gap 1e-12", f"{tag}: tail gap")
    match = re.search(r"precheck_clower (\[[^\n]+?\])", graph_log)
    require(match is not None and lower_endpoint(match.group(1)) > 0,
            f"{tag}: graph tail lower bound")

manifest = ROOT / "SHA256SUMS.txt"
expected_manifest_paths = {
    "PARENT-SHA256.txt",
    "README.md",
    "results/even-checker.log",
    "results/even-contract.log",
    "results/even-family.log",
    "results/even-graph.log",
    "results/even-safe.log",
    "results/log11_even_K400_N500_clean_sharp_family_500_850_full393_q8000_arb_256.npz",
    "results/log11_odd_K400_N500_clean_sharp_family_500_850_full394_q8000_arb_256.npz",
    "results/odd-checker.log",
    "results/odd-contract.log",
    "results/odd-family.log",
    "results/odd-graph.log",
    "results/odd-safe.log",
    "sources/export_log10_even_yukawa_family_parseval_data.py",
    "sources/extract_log10_even_K400_full_contract_base.py",
    "sources/generate_log10_endpoint_contract_cache.py",
    "sources/log10_even_K400_sharp_safe_arb.py",
    "sources/log8_block_certificate_body.py",
    "sources/verify_first_mixed_band_endpoint_log5.py",
    "verify_frozen_results.py",
}
lines = manifest.read_text().splitlines()
require(len(lines) == 21, "manifest entry count")
entries = []
for line_number, line in enumerate(lines, 1):
    fields = line.split(maxsplit=1)
    require(len(fields) == 2, f"malformed manifest line {line_number}")
    digest, rel = fields
    rel = rel.removeprefix("./")
    entries.append((digest, rel))
manifest_paths = [rel for _, rel in entries]
require(len(set(manifest_paths)) == len(manifest_paths),
        "duplicate manifest path")
require(set(manifest_paths) == expected_manifest_paths,
        "manifest path set")
for digest, rel in entries:
    target = ROOT / rel
    require(target.is_file(), f"missing manifest target {target}")
    require(hashlib.sha256(target.read_bytes()).hexdigest() == digest,
            f"hash mismatch {target}")

print("log11 round-1 frozen results: metadata, pivot counts, signs, hashes passed")
