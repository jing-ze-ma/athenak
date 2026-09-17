---
name: general-eos-table-linear-wave
description: "How general_eos=table is regression tested (linear-wave convergence), and the host-vs-device constraint that shaped the pgen changes"
metadata: 
  node_type: memory
  type: project
  originSessionId: 2cb708dd-0a5f-46c0-b912-dd406f996581
  modified: 2026-08-13T12:58:58.552Z
---

Commit `4a59a1c9` on `general-eos` (2026-08-13) closed items 1 and 2 of
[[general-eos-project]]'s resume list.

**The constraint that drives everything here: the tabulated EOS CANNOT be evaluated from
host code.** `EOSTable::tbl` is a `DvceArray3D`, so on a GPU build any host-side
`eos.Pressure/Gamma1/EnergyFromPressure` reads device memory. It only "works" on a CPU
build by accident. So every EOS call in a pgen must sit inside a `par_for`. When a pgen
genuinely needs ONE number on the host, use the new
`pgen_eos::HostGamma1FromP(eos,d,p)` in `src/pgen/pgen_eos_utils.hpp`, which launches a
one-element kernel and `deep_copy`s the answer back. Do not add host-side EOS calls.

**`shock_tube` / `linear_wave`:** `<problem>/pl,pr,pgas` are now real pressures under
`general_eos=table`; the p->e conversion moved into the kernels. The ideal branches keep
`p/(gamma-1)` verbatim (NOT `p*igm1` — they differ in the last bit, verified), and hydro +
MHD linear waves and a cooling shock tube were confirmed bitwise unchanged by rebuilding
the stashed tree and comparing `%24.16e` `.tab` output.

**`linear_wave`'s eigensystems now take `gam1`.** In the PRIMITIVE variables
`(d,vx,vy,vz,P[,By,Bz])` the ratio of specific heats enters the hydro and MHD eigenvectors
ONLY through the sound speed — in MHD the remaining `gm1`s sit in
`bt_starsq = (gm1-(gm1-1)y)btsq` and the `bet*_star`, which collapse to 1 at the `y = 1`
this generator passes. So substituting `Gamma_1` for `gamma` is exact, not an
approximation. (`h`/`hp` are computed but unused.)

**The new tests** (`tst/scripts/{hydro,mhd}/*_general_eos_table.py`, inputs
`inputs/tests/linear_wave_{hydro,mhd}_geneos_table.athinput`) are CONVERGENCE tests, since
there is nothing to compare a real EOS against. cgs background: rho = 1e-6 g/cm^3,
p = 6.5e5 erg/cm^3 (~1e4 K), in the H partial-ionization zone where **Gamma_1 = 1.30**,
not 5/3. Measured order 1.91-2.00 over nx1 = 64/128/256. Each also runs the same input in
`general_eos=gamma` mode and requires the error levels to differ (~2x), so the test cannot
pass on a silent fallback.

**Trap that cost a debugging cycle:** the linear-wave eigenvectors are NOT normalized —
the density component is 1, so `<problem>/amp` is an ABSOLUTE density perturbation. The
usual code-unit tests hide this with `dens = 1`. In cgs with `dens = 1e-6` you need
`amp = 1e-10` for a 1e-4 relative wave; `amp = 1e-4` gives velocities of 1e7 cm/s and a
dt of 1e-10 s.

Wave speed is a free check on the EOS: the pgen sets `tlim` to one wave period, so
`cs = L/tlim`. Table hydro gave 1088.09 s (cs = 9.19e5, Gamma_1 = 1.30) vs gamma mode
960.77 s (cs = 1.04e6, gamma = 5/3); the MHD fast speed matched the analytic
`cf` to all printed digits.
