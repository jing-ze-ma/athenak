---
name: ck-mhd-b3-run
description: ck_mhd_b3 -- MHD version of the current dhj production setup with EOS resistivity at bbot = 3 G. Built and SMOKE-TESTED, not launched
metadata:
  type: project
---

`/viper/u2/jinma/ATHENAK/bench/ck_mhd_b3`, set up 2026-08-27 at the user's request
("test mhd version of the problem with eos resistivity, maybe just choose bbot = 3 first...
but with updated point mass gravity and all the updates").

It is [[ck-grav-prod-run]]'s input -- point-mass gravity, correlated-k RT, tabulated general
EOS with metal ionization, the polynomial radial stretch, the measured floors -- with **six
lines changed and nothing else**, so hydro-vs-MHD is a clean comparison at fixed everything:
`<hydro>` -> `<mhd>`, `rsolver` hllc -> hlld, `ohmic_resistivity = eos`, `max_eta = 1.0e13`,
`bbot` 1.0e1 -> 3.0e0, output3 `hydro_w` -> `mhd_w_bcc`.

**No rebuild was needed.** AthenaK picks the physics module from which input blocks exist,
and `deep_hot_jupiter_rt` already carries `pmhd` branches in the IC, the hydrostatic outer
BC, the source terms and all three RT schemes. Same `athena` snapshot, HEAD `1e19d4e7`.

**`ohmic_resistivity = eos` needed no extra configuration**, unlike [[xe-resistivity-long-runs]]
which had to bolt a `<units>` block and the `eos_*` composition keys onto an ideal gas.
`resistivity.cpp` fatals unless the EOS carries an electron fraction, and
`eos = general` + `general_eos = table` + `eos_metal_ionization = true` is already exactly
that. `max_eta = 1e13` and `use_rkg_sts = false` are carried from the xe_long measurements.

## SMOKE 11110445 PASSED -- and the field costs nothing in dt

Against the hydro smoke 11105755 at the identical t = 3217.5 s: **dt at cycle 0 is 10.8648 s
in BOTH**, mass drift +3.118e-3 vs +3.128e-3, energy +3.556e-3 vs +3.503e-3. So at 3 G the
Alfven speed never approaches binding, and the drift is the pre-existing initial transient,
not the field. This confirms by measurement the one thing I had only estimated: `dfloor` here
is 145x lower than xe_long's (5e-14 vs 7.26e-12), which makes v_A at the domain top ~8 km/s
against dr 2.18e8 cm, i.e. a ~270 s Alfven CFL against the ~11 s hydro dt.

**Caveat, do not misreport:** the floor/C2P counters were NOT sampled -- `<output2>` has
dt = 3.05e4 and the smoke only reached t = 3217, so `dhj.log` has only its header. "No floors
fired" is NOT established.

**NOT LAUNCHED.** `sbatch submit_viper.sh` from that directory is the remaining step.
