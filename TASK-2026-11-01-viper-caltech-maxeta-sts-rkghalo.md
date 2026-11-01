# TASK 2026-11-01 viper -> Caltech (USER): (A) max_eta test + STS vs explicit; (B) gate of the exact RKG halo saving (rkghalo-test-1101)

The USER asked for both, using your idle GPUs.

## A. Does the max_eta cap set the cost, and can it go lower? And should we use STS (RKG) or plain explicit?

### Why
- On the 10 G state the binding resistive dt is set by cells at the cap, max_eta = 1e15 cm^2/s.
- With Saha x_e, eta >= 1e15 in 9-12 % of cells at 1e-4..1e-1 bar (night 1e-2..1e-1 bar and the terminators; the x_e tracer does not change these cells).
- The global min h^2/eta is 0.246 s Saha / 0.253 s tracer, i.e. h ~ 160 km at eta = 1e15. RKG stages scale as sqrt(max_eta).
- In capped cells Rm << 1 whatever the exact eta, so B slips through either way. With uL ~ 1 km/s x 300 km = 3e12 cm^2/s:
  - cap 1e14 -> Rm(cap) 0.03, ~3x fewer stages;
  - cap 1e13 -> Rm(cap) 0.3, ~10x fewer stages (borderline).
- Explicit resistivity: at cap 1e15 the resistive dt (~0.25 s) is ~15-20x below the hydro dt, so RKG wins. At 1e13 the resistive dt (~25 s) is ABOVE the hydro dt (~4-5 s), so explicit needs no super-stepping at all and has no splitting error. 1e14 is in between.

### Arms
- Binary: fin51-r-fix head (65658289c or newer).
- Deck: the r1m deck with its production keys (fac 4, outflow1, cs_seam_flux shared, am_hst_split false).
- Start: one developed 10 G state, e.g. your p_f4 0.5-rot rst or newer. All arms start from the same rst.
- Duration: 1.5 rot each unless noted.

| arm | max_eta | resistive integrator |
|---|---|---|
| R15 | 1e15 | RKG as in production (rkg_dt_safety 0.5, rkg_lean): the reference |
| R14 | 1e14 | RKG, same keys |
| R13 | 1e13 | RKG, same keys |
| X14 | 1e14 | explicit (use_rkg_sts = false) |
| X13 | 1e13 | explicit |
| X15 | 1e15 | explicit, ~300 cycles only, for cost |
| C (optional) | 1e15 | R15 with a different rank count or layout, as the chaos-spread reference |

### Report
- Cost: s/cycle, mean dt, RKG stages/cycle, s/rot, for every arm.
- At 0.5 / 1 / 1.5 rot, against R15 (and C):
  - horizontal-mean and equatorial u_phi at 1e-4, 1e-3, 1e-2, 1e-1 bar;
  - day / night / terminator T at the same levels;
  - B_rms per band and E_B;
  - Joule heating or ohm_drag if your deck outputs it;
  - KE, total-E drift;
  - event counts.

### Decision rule (we will decide with the USER)
- A lower cap is acceptable if its differences from R15 are within the R15-C spread, or < 2 % in winds and T at p <= 0.1 bar, with no systematic trend in B_rms at the Fe level (1e-5..1e-3 bar).
- Between RKG and explicit at the same cap: explicit wins if it is cheaper, or within ~20 % of the cost. It has no splitting error.

## B. Exact RKG halo saving: fork/rkghalo-test-1101 @ 0594b145d (viper apudev queue blocked)

### What it does
- New key `<mhd>/rkg_halo_ien_only` (default false). In RKG super-stages 2..s it exchanges only u0(IEN); every other variable is unchanged during RKG.
- An exact fallback: if ConToPrim changed any non-IEN value in an active cell or an exchanged ghost on any rank (MPI_Allreduce), the next stage does the full exchange.
- It does nothing with rkg_lean (FATAL if combined) or without use_rkg_sts. So it is for NON-lean runs only. Your production uses rkg_lean, so this mainly matters for the explicit/RKG choice above and for any bitwise-required run.
- Expected: ~43 of 263 ms/cycle on viper (s10, non-lean), if the fallback is rare.

### Tests (build 0594b145d as usual)
1. **Smoke:** an r1m restart, ~10 cycles, with the production keys but rkg_lean = false, plus `mhd/rkg_halo_ien_only=true`. No FATAL.
2. **Bitwise:** same binary, same restart, +600 s (or ~200 cycles), key off vs on. Write hst every cycle (output1/dcycle=1) and a final rst.
   - rst bytes after `<par_end>`: identical.
   - dhj.mhd.hst: identical.
   - dhj.user.hst: within 1e-12 (the Mfloor/Efloor columns are GPU-atomic noise).
3. **Fallback count:** grep `rkg_halo_ien_only:` in the key-on log. Report N together with the number of RKG stages and cycles.
4. **Cost:** key off/on interleaved, 3 repeats, >= 300 cycles each, s/cycle.

Please report A and B in a NOTE. Viper keeps its own B jobs queued and will cross-check.
