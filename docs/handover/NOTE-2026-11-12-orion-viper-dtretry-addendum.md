# NOTE 2026-11-12 Orion -> viper: dtretry addendum (34d198efa, m1.vgd_hr in the minimal list)

Build 34d198efa (PROBLEM=he_star_m1, Release, MPI, -ffp-contract=off): athena md5 c9a5f38d2697995838a33a07a3ad197e.
run2.sh unchanged, RUNS=runs2, same tables as before; host orion003, 10 arms x 3 ranks.

**Result: the fix works on CPU.** CN vs CF and CTN vs CTF are now **IDENTICAL** (00002-00004); they DIFFERED at 3ec401877.
All 10 arms rc=0, 0 FATAL. V1..V6 are all IDENTICAL to CN: with vgd_hr already in the list, no extra group is needed.

## gate2.out (whole)
```
CN rc=0 fatal=0 
CTN rc=0 fatal=0 
CTF rc=0 fatal=0 
V4 rc=0 fatal=0 
V1 rc=0 fatal=0 
V3 rc=0 fatal=0 
V6 rc=0 fatal=0 
CF rc=0 fatal=0 
V5 rc=0 fatal=0 
V2 rc=0 fatal=0 
== CN vs CF
hegiant.00002.rst IDENTICAL
hegiant.00003.rst IDENTICAL
hegiant.00004.rst IDENTICAL
== CTN vs CTF
hegiant.00002.rst IDENTICAL
hegiant.00003.rst IDENTICAL
hegiant.00004.rst IDENTICAL
== CN vs V1
hegiant.00002.rst IDENTICAL
hegiant.00003.rst IDENTICAL
hegiant.00004.rst IDENTICAL
== CN vs V2
hegiant.00002.rst IDENTICAL
hegiant.00003.rst IDENTICAL
hegiant.00004.rst IDENTICAL
== CN vs V3
hegiant.00002.rst IDENTICAL
hegiant.00003.rst IDENTICAL
hegiant.00004.rst IDENTICAL
== CN vs V4
hegiant.00002.rst IDENTICAL
hegiant.00003.rst IDENTICAL
hegiant.00004.rst IDENTICAL
== CN vs V5
hegiant.00002.rst IDENTICAL
hegiant.00003.rst IDENTICAL
hegiant.00004.rst IDENTICAL
== CN vs V6
hegiant.00002.rst IDENTICAL
hegiant.00003.rst IDENTICAL
hegiant.00004.rst IDENTICAL
wall 408 s
```

Restart trees: runs (3ec401877 gate, ~30 GB) and runs2 kept on Orion under /orion/ptmp/jinma/dtretry_1112; I will delete them on your OK.
