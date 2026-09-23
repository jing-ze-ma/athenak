# Is a deep-hot-Jupiter production restart bitwise? (09-23)

Run tree: /viper/ptmp2/jinma/rstjolt_0923 (CPU, MPI, 6 ranks, login node nice 10).
Binary: rt-integration be02c647, -D PROBLEM=deep_hot_jupiter_rt, Athena_ENABLE_MPI=ON,
Release, gcc/14 openmpi/5.0 (build_cpu in the wt_rstjolt worktree).
Inputs: the production inputs with hst every cycle at %.17e (`dcycle = 1`,
`data_format = %.17e`) and bin hydro_u every cycle (hydro).
Restarts, read in place: hydro bench/cs_hyd4_prod/rst/dhj.00131.rst (t = 1.99777e7 s,
rotation 65.5, cycle 1019222); MHD bench/cs_mhd_prod4/rst/dhj.00100.rst (t = 1.525e7 s,
rotation 50.0, cycle 975801).

Test: A = 20 cycles straight; B1 = 10 cycles, rst at the end; B2 = restart from B1's rst,
10 more cycles. Compare A and B2: hst rows, bin, rst payload after `<par_end>`
(`rstcmp.py`). A run twice is bitwise (determinism control).

## Before the fix (be02c647)

| | hst (20 rows, %.17e) | rst payload | bin |
|---|---|---|---|
| hydro | bitwise | BITWISE | bitwise |
| MHD | differs from the 1st cycle after the restart (1-mom 5e-14 rel) | DIFFER | - |

MHD: the rst B1 wrote is bitwise equal to A's rst at the same cycle; one cycle later
(`rstdiff.py`) everything above i_act ~ 37 differs. The in-memory state at the start of
the first restarted cycle (a temporary dump hook in Driver::Execute, not committed; `dcmp.py`)
differs ONLY in w0(IEN), by <= 7e-14 relative in ~1e-2 of the cells (i >= 41); u0, b0,
bcc0, wtemp, wder, eta_b are identical (u1, b1, efld, wbq0 are per-stage scratch). Hydro
has the same w0(IEN) difference (<= 5e-14 in ~600 cells) but it did not reach u0 within
10 cycles.

Missing state: the internal energy of the last inversion. With etotgrav, ConToPrim does
RemoveGravEtot -> inversion -> AddGravEtot on u0 (hydro/hydro_tasks.cpp Hydro::ConToPrim,
mhd/mhd_tasks.cpp MHD::ConToPrim). The file stores E' = (E - rho phi) + rho phi; the restart's
first (frozen) inversion uses E' - rho phi, which is not the E - rho phi the straight
run used, so w0(IEN) = E' - rho phi - KE - ME differs in the last bits
(coordinates/gnomonic_raisevel_frozen.cpp GnomonicEquiangleRaiseVelMHDFrozen, w0(IEN) line).

Size in production terms (MHD, T = wtemp, `jolt.py`): 1 cycle after the restart,
max |dT|/T = 3.1e-4 in the top active cell (i_act = 127); 117 of 786432 cells above 1e-6;
below the top cell <= 2e-9 (i_act 56-120) and 0 below i_act ~ 40. After 10 cycles:
max 3.1e-5 (top cell), 148 cells > 1e-6, <= 6e-9 elsewhere. This is roundoff seeding
amplified in the top cell, not a jolt of the 5 K / 1.3e-3 kind: the wellposed A1 jolt
(tests_ck_implicit/wellposed/README.md) does not appear with the production switches and
comes from something else.

## Fix (branch rst-jolt)

A marked restart block `EINTRST1` (kEintRstMagic, pgen/pgen.hpp) behind the other marked
headers, with one slab per general-EOS module (hydro, then MHD) of w0(IEN), as the last
slabs of each MeshBlock record (outputs/restart.cpp). The reader (pgen/pgen.cpp) loads it
into Hydro/MHD::eint_rst only when the full general-EOS tail is there (the frozen first
conversion). That conversion then puts eint_rst into w0(IEN) and puts back u0(IEN) as the
file had it, before SetResistivity. Files without the block still read, with a warning
(the old behaviour). The continuous run is unchanged: hst of the straight 20-cycle runs is
byte-identical to be02c647, and so is the old-layout part of every rst record.

After the fix (A vs B2): MHD rst payload BITWISE, 11 hst rows byte-identical; hydro rst
payload BITWISE, 11 hst rows byte-identical.

Note: an older binary (the running productions) cannot read a file with the new block.
The productions write and read their own format, so this does not affect them.

Scripts: run.sh (one run), cmp.py (hst by time), rstcmp.py (rst payload), rstdiff.py
(per-variable rst diff, MHD layout), jolt.py (|dT|/T from wtemp), dcmp.py (in-memory dumps).
