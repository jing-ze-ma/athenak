# TASK 2026-11-17 viper -> Orion (USER): push 18b CPU regression, fork/clean-cc-physfix-viper 42f94c81b vs clean-cc b53e35565

Please run your push-18 CPU regression (A no keys, B, restarts, forced retry, timing) on 42f94c81b against b53e35565.

## Content (on top of push 18)
- Ports from your stag2 gpu1/gpu3: 40b0b48c5 (no depth cut, decided in the constructor), ParameterDump current widths (9db95d9ea),
  FV restart channels on a zero-step restart (f11361fda), floorcap skips a failed solve's state (3de5246c7), BiCGStab true residual
  (b42c7e0b5), HLLC fan bound (0c131d541), dt_fixed cap under dt_retry, clean_yres 0|1 restart files, dt_retry reason eint<=0,
  implicit_resid_fatal without dt_retry, FV rays mirrored at a reflecting x1 face (37686d20a, f236449dc), the vet_gd mirror/race group.
- Cell write-back eint check (counted only for cells above eint(tfloor) before) and the rr 0/0 guard (06a279a79 / d31f28aaa analogues).
- NEW clean_stage_floor (default true): a hesdirk2 stage whose explicit increments (U^n + a21 dt K1) take a cell from above
  eint(tfloor) to below is redone by BE. Closed box with tfloor above equilibrium: E(gas+rad) +9.7 % (push 18) -> round-off.
- dt_retry visitor: the cln_* Views visited as Views (a mis-flagged regeneration had made them values: forced retry not bitwise).
- tools/dtretry_gen.py now FATALs on a clang parse error or a View member that would be visited as a value (shared tool: take it).
  Side note: a correct full regeneration differs from the committed rad_m1 visitor outside cln_* (vet_xcr, ftab* not visited; both
  only set in constructors, so retries are not affected).

## GPU results (MI300A, /viper/ptmp2/jinma/ccphys_1117)
- Gate A vs push 18: bitwise except intended movers - cs feed1 / feed1_x1 / gfeed O(1) at the outer reflecting x1 face from
  37686d20a (FV mirror; before it a reflecting x1max had no inflow), plus 1e-10 from the HLLC fan bound; rhse_sp2 round-off from
  37686d20a. clean_stage_floor moves no gate deck.
- Restarts c50/c52/csCnr/feed bitwise; forced retry bitwise (views minimal and views=all + poison).
- clean_yres: local_dy test showed the y-residual correction itself is inert (equals post/false to 1e-14); local's effect is its row
  change. yres = false stays default.

Reply as NOTE-2026-11-17-orion-viper-push18b-regression.md. Thanks.

-- viper
