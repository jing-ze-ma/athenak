"""nbutil.py: the two helpers the notebooks share (Python standard library + numpy).

suite(test, exe, work): run one test of the suite (run_tests.py --only TEST --keep) unless
  WORK/TEST.json already exists, and return its measured values (the --record JSON).
verdict(test, res): print the suite's PASS / FAIL line with the criterion of
  ref/references.json; returns True on PASS.
"""
import json
import os
import subprocess
import sys

SUITE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, SUITE)
import vett  # noqa: E402,F401  (re-exported for the notebooks)
import run_tests as rt  # noqa: E402


def suite(test, exe, work, threads=4):
    rec = os.path.join(work, test + '.json')
    if os.path.exists(rec):
        print('loading the existing run', os.path.relpath(os.path.join(work, test)))
    else:
        cmd = [sys.executable, os.path.join(SUITE, 'run_tests.py'), '--exe', exe,
               '--only', test, '--keep', '--workdir', work, '--threads', str(threads),
               '--record', rec]
        print('running:', ' '.join(cmd), flush=True)
        r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           universal_newlines=True)
        print(r.stdout)
        if not os.path.exists(rec):
            raise RuntimeError('the suite wrote no record (see the output above)')
    return json.load(open(rec))[test]


def verdict(test, res):
    bad = rt.check(test, res, rt.REF[test]) if res.keys() - {'wall'} else ['the run failed']
    crit = {k: v for k, v in rt.REF[test].items() if not k.startswith('_')}
    print('%s %s   criterion %s' % ('PASS' if not bad else 'FAIL', test, json.dumps(crit)))
    for b in bad:
        print('   <- ' + b)
    return not bad
