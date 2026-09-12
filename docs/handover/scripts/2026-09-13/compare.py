#!/usr/bin/env python3
"""Compare two run_cases.py result files case by case.

Prints, for every case present in both: the PASS/FAIL verdict on each side and, for
each metric the driver recorded, the reference value, the new value and the relative
difference.  Last-bit drift between a CPU and a HIP build is expected, so the verdict
(which applies the tests' own thresholds) is the gate; the numbers are here to show
HOW far apart the two builds are.
"""

import re
import sys

_METRIC = re.compile(r"([A-Za-z_][A-Za-z_0-9]*)=(\(?[-+0-9eE.,naif]+\)?)")


def parse(path):
    """label -> (verdict, {metric: [values]})"""
    out = {}
    with open(path) as fp:
        for line in fp:
            if not (line.startswith("PASS ") or line.startswith("FAIL ")
                    or line.startswith("SKIP ")):
                continue
            verdict, rest = line.split(None, 1)
            label = rest.split()[0]
            metrics = {}
            for name, val in _METRIC.findall(rest):
                if name in ("t",):
                    continue
                metrics.setdefault(name, []).append(val)
            out[label] = (verdict, metrics)
    return out


def relnum(a, b):
    try:
        fa, fb = float(a), float(b)
    except ValueError:
        return None
    if fa == fb:
        return 0.0
    denom = max(abs(fa), abs(fb), 1e-300)
    return abs(fa - fb) / denom


def main():
    ref_path, new_path = sys.argv[1], sys.argv[2]
    ref, new = parse(ref_path), parse(new_path)
    print(f"# reference : {ref_path}   ({len(ref)} cases)")
    print(f"# new       : {new_path}   ({len(new_path and new)} cases)")
    only_ref = sorted(set(ref) - set(new))
    only_new = sorted(set(new) - set(ref))
    common = [k for k in new if k in ref]

    agree = differ = 0
    for label in common:
        rv, rm = ref[label]
        nv, nm = new[label]
        marks = []
        for name, vals in nm.items():
            rvals = rm.get(name)
            if not rvals:
                continue
            for i, v in enumerate(vals):
                if i >= len(rvals):
                    break
                d = relnum(rvals[i], v)
                if d is None or d == 0.0:
                    continue
                marks.append(f"{name}[{i}] {rvals[i]} -> {v} (rel {d:.2e})")
        same = (rv == nv)
        if same and not marks:
            agree += 1
            continue
        differ += 1
        print(f"{'VERDICT ' + rv + '->' + nv if not same else 'metrics differ':<24}"
              f" {label}")
        for m in marks[:8]:
            print(f"      {m}")

    print()
    print(f"# {len(common)} cases in common: {agree} identical verdict and numbers, "
          f"{differ} with a differing verdict or metric")
    for label in only_new:
        print(f"# only in {new_path}: {label} ({new[label][0]})")
    for label in only_ref:
        print(f"# only in {ref_path}: {label} ({ref[label][0]})")
    nfail = sum(1 for k in new if new[k][0] == "FAIL")
    print(f"# {new_path}: {nfail} FAIL of {len(new)}")


if __name__ == "__main__":
    main()
