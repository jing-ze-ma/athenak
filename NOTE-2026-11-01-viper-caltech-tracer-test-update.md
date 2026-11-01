# NOTE 2026-11-01 viper -> Caltech (USER): tracer tests, build 89213fcb3 instead of 2b72d789; adjusted comparison rules

Update to TASK-...-tracer-gpu-tests (bf54375e6).

## 1. Build fork/tracers-test-1101 @ 89213fcb3
- On 2b72d789 the new hst columns Mdrain/Edrain were always ~0: a local `ke` shadowed the loop bound.
- The drain itself and its booking into Mfloor/Efloor were already correct. Only those two diagnostic columns change.
- If T0 already ran on 2b72d789, its key-off result still holds: the fix touches only key-on diagnostics.

## 2. T0 key-off rule on GPU
- **Compare exactly:** the rst data and dhj.mhd.hst.
- **dhj.user.hst:** compare with the floor columns (Mfloor/Efloor/Efloor_rt) at <= 1e-12 rel, or excluded. Those columns are summed with GPU atomic_add (c2p_track), so the summation order varies from run to run.
- **Viper evidence:**
  - b935e145 vs the branch, s10 restart +600 s, plain and with r1m keys: rst and mhd.hst identical; user.hst differs only in those three columns, at 6e-16..4e-13.
  - On CPU, everything is identical.
  - A b935-vs-b935 GPU repeat is queued to confirm that this is base behaviour.

## 3. T2 restart continuation: how to compare
- Use cycle-matched runs (time/nlim), not tlim. The straight run goes to cycle N; the run restarted from its mid rst goes to the same N.
- Compare the DATA tail of the two final rst files, e.g. `cmp <(tail -c $((size-20000)) a.rst) <(tail -c $((size-20000)) b.rst)`. Whole-file cmp fails even when the state is identical, because:
  - the header carries the input text;
  - the implicit-ck restart section changes size with the ck cadence phase.
- On viper CPU (both tracers), 150 + 150 cycles vs 300 cycles gives an identical data tail.

## 4. Add to T1
- Mdrain [g/s] and Edrain [erg/s] are the last two flux-hst columns. Gross floor = Mfloor (net) + Mdrain.
- **Expected:** at the default margin 2 / fmax 0.01, net Mfloor falls well below gross, and the mass ledger closes as in the off arm.
- **Viper CPU stress test** at margin 1 / fmax 0.1: the ledger closes to 0.2 % of Mdrain. It also shows the pump that margin 2 prevents: the floor re-adds what the drain removes (gross ~= Mdrain).
- **Note:** a late-time state (s10, Mfloor 1.3e15) drains little in 900 s, because its floor gas sits at rho < 1.5 rho_fl. The from-t=0 r1m arms are the meaningful test.

## 5. Viper GPU key-on smoke
s10 restart without tracers, r1m keys, +900 s. Arms off / xe / fl / both: all clean. Restarts both->fl and fl->both: clean. Mass-ledger residual ~1e10 g/s (hst sampling).
