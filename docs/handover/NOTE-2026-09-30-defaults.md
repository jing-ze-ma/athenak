# NOTE 2026-09-30: four code defaults flipped (user decision 09-30)

Branch `defaults-0930` (commit b7f0d8a1), merged into rt-integration. Source: the class-A entries and the
`sponge_bottom` flag of `AUDIT-2026-09-30-default-candidates.md` (this directory).

| key | old default | new default | where |
| --- | --- | --- | --- |
| `problem/ck_impl_maxit` | 8 | 24 | src/pgen/deep_hot_jupiter_rt.cpp:739 (GetOrAddInteger); src/utils/two_stream_column_ck.hpp:114 |
| `problem/ck_impl_conserve` | 0 | 1 | src/pgen/deep_hot_jupiter_rt.cpp:926-928 (read only if named); src/utils/two_stream_column_ck.hpp:181 |
| `problem/bbot` | none (FATAL if absent) | 0.0 | src/pgen/pgen.cpp:92 (ReadHotJupiterBbot, GetOrAddReal); src/pgen/deep_hot_jupiter_rt.cpp:2282 |
| `problem/sponge_bottom` | true | false | src/pgen/pgen.cpp:132 and :212 (both ProblemGenerator ctors); src/pgen/pgen.hpp:74 |

`bbot` + `bbot_gauss` together is still FATAL; `bbot_gauss` unchanged. The dhj pgen now prints one startup line
`deep_hot_jupiter_rt: ck_impl_maxit N, ck_impl_conserve N, sponge_top N, sponge_bottom N`.

## Restarts

`ParameterInput::GetOrAdd*` returns the stored value when the key exists; a restart loads the parameter dump
embedded in the rst, which holds every key read with GetOrAdd in the original run. So `ck_impl_maxit`,
`sponge_bottom` and `bbot` keep the values of the old run on restart. `ck_impl_conserve` is read only when named
(it is never added to the rst), so an old restart that did not name it switches to 1: reporting only (the
`Lir_top`/`Etot_top` user-hst columns), the solution is bitwise.  A key that the restart does not carry cannot be
set on the command line (`-r ... problem/ck_impl_conserve=0` is FATAL "not found").

## Gates (viper apudev, ROCm 7.2, tip = cba4e797 binaries = fork/rt-integration baed412e src)

- (a) w121prod_0929/w1x/rst/dhj.00600.rst (rot 300.0, 1x hydro), 100 cycles, 2 ranks (job 12037161):
  rst and bin BITWISE tip vs new, hydro.hst identical; user.hst differs in Lir_top, Etot_top (conserve 0 -> 1,
  that restart does not carry the key) and in the last digits of Efloor/Mfloor/Efloor_rt, which differ equally
  between two tip runs (job 12037467, run-to-run reduction noise).
- (b) fresh start from w121prod_1x.athinput with the four keys removed (job 12037161): runs 20 cycles, log
  `ck_impl_maxit 24, ck_impl_conserve 1, sponge_top 1, sponge_bottom 0`, `bbot = 0`; the rst records
  bbot 0, sponge_bottom 0, ck_impl_maxit 24.
- (b') same input with the OLD values named (maxit 8, conserve 0, bbot 0.0, sponge_bottom true), tip vs new,
  20 cycles (job 12037546): rst, bin, hydro.hst BITWISE; user.hst only Efloor/Efloor_rt last digits (noise).
- (c) He box (hebox_cfl2_0927 input, C15 rst t = 3800, cfl 0.9) 30 cycles, box_convection tip vs new: rst,
  bin, both hst BITWISE (box_convection reads none of the four keys; pgen.cpp reads them only if hot_jupiter).
- (d) tst: test_rad_dhj_ck_cpu, test_rad_dhj_srclim_cpu, test_rad_dhj_ck_mpicpu PASSED (none has reference
  numbers that depend on these keys; no test input changed). cpplint: no new findings in the four files.

Gate tree: /viper/ptmp2/jinma/defaults_0930/gate (gate*.sub, cmp.py).

## Restoring the old behaviour

Name the old values in `<problem>`: `ck_impl_maxit = 8`, `ck_impl_conserve = 0`, `sponge_bottom = true`
(`bbot` needs nothing: an input that names it behaves as before). Gate (b') shows this is bitwise to the old code.
