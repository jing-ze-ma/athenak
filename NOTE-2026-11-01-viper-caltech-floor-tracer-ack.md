# NOTE 2026-11-01 viper -> Caltech (USER): floor tracer taken over; x_e tracer status; key recommendation for r1m (preliminary)

Reply to NOTE-...-noinfall-fac-smokes (a2830c93a) and TASK-...-floor-tracer-to-viper (b42671501).

## Floor-mass tracer: taken over (started now)
- It will be built next to the x_e tracer, on one branch with one restart path. Restarts will work from a file with none, one or both scalars (missing x_e -> Saha, missing floor -> 0).
- Your WIP 45010a2b is the starting point. We will report what we change.
- Gates: your list (key off bitwise CPU C8x112 + GPU; key on clean, restart without the scalar, 0 <= s <= rho), plus an explicit mass/energy ledger and a viper GPU smoke on the r1m keys (fac 4, outflow1, rkg_dt_safety 0.5, rkg_lean), compared with your p_f4 baseline.

## x_e tracer: built, final validation running (local branch xetracer-1101, not yet on fin51-r-fix)
- **Key:** `<mhd>/xe_tracer`, one MHD passive scalar z = x_e.
- **Source:** dn_e/dt = alpha(T)(n_eq^2 - n_e^2), Na+ radiative recombination, backward Euler in closed form, once per cycle. z replaces the Saha x_e in eta_O and in the ambipolar rho_i.
- **Gates so far:**
  - Key off: GPU restart files identical on s10 (RKG) and b100 (explicit).
  - CPU unit tests: relaxation and advection.
  - Restart from a file without the scalar starts from Saha exactly.
  - The single-carrier model matches our full offline chemistry to a median 1.0-1.2x at 1e-6..1e-1 bar.
- **Cost:** +6.0 %/cycle on MI300A, s10 restart. A kernel profile is running to see what can be saved.
- **Running:** a 2-rotation s10 run with the tracer, to compare night F against the offline 80-800x and to measure the change in night eta_O, Rm and the resistive dt. Then rebase onto fin51-r-fix (b935e145 + ambipolar) and NOTE.

## Ambipolar: being merged into fin51-r-fix
`<mhd>/ambipolar`, default off = bitwise. The merge is being re-gated on b935e145 (rkg_lean / rkg_dt_safety / noinfall).

## Recommended r1m keys (preliminary; final in the tracer NOTE)
- **xe_tracer = true** once it is on fin51-r-fix.
- **ambipolar: off for r1m.**
  - It is negligible on the day side and terminator (AD/advection < 1e-3).
  - At night it was ~1 at 1e-6..2e-4 bar only with Saha x_e. eta_A scales as 1/rho_i, so with the tracer it drops by the same F (80-800x).
  - On a restart it shortened the explicit dt (0.82 -> 0.63 s).
- **valfven_noinfall_fac = 4,** as you recommend, plus floor_tracer when ready.
