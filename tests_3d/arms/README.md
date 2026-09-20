# he4_presn on the cubed sphere: the 3-D SMOKE grid, arms A and B
(2026-09-17, viper, branch `he4-presn-global`, binary `../../build_gpu_rg/src/athena`)

The 1-D gate of `../README.md` is red and cannot be made green for this star: a
super-Eddington, radiation-dominated FeCZ has no hydrostatic radiative-equilibrium state,
so 15-40 % of `L` *must* be carried convectively.  Per the brief the decision was to go
straight to 3-D from the MLT-built ic
(`/viper/u2/jinma/ATHENAK/bench/hestar_presn/ic_he4_presn_sph.txt`) and gate on the 3-D
behaviour instead.

**Both inner-boundary arms FAIL.  Arm A (open) dies at t = 549 s = 0.117 turnover,
arm B (wall + base flux) at t = 1723 s = 0.366 turnover.  Arm C was not run: its
precondition (a survivor past one turnover) was not met.**

## Setup

`inputs/hydro/he4_presn_cs.athinput` as committed: nx1 = 96 stretched from 0.50 R,
nx2 = nx3 = 32 per panel, `meshblock 96x16x16` = 24 blocks, 2 MI300A (`srun -n 2`,
`--gres=gpu:2`), `tlim = 1.4115e4` s = 3 turnovers (turnover = 4.705e3 s), hst and
`rt_profile` every 47 s, `rt_surface` every 94 s, `bin` every 0.5 turnover, `rst` only at
exit.  All production switches on: mode 3 (`pcr`, `mixed 2`, `warm 1`), `rt_rad_force`,
`rt_strang`, the implicit ADI transverse operator, FOFC, both sponges, `wb_cache_every = 0`.
Seed: `vpert = 1e-3` c_s on `v1`, tau window 20..300, `vpert_lmin/lmax = 2.5e11/7.4e11`
(l ~ 2-6, about 2 Hp_FeCZ).

One input change was needed and is committed with this note: an `<output5>`
`file_type = log` block, so the run writes the **event log** (`he4.log`) the gate reads.
Nothing else in the input or in `src/` was touched -- no bug was found that called for it.

| arm | inner boundary | `rad_flux_inner` | `rt_bottom_flux` |
| --- | --- | --- | --- |
| A | `open` (deep-adiabat relaxation) | 0 | false |
| B | `wall` | 1.305278e15 = L/(4 pi r_in^2) | true |

`chain.sh` is the submit script (`--export=ALL,ARM=<A|B>,WALL=hh:mm:ss`); it restarts from
`rst/` if one is there, so an arm can be chained across short jobs.  `gate.py <armdir>`
builds the tables below from `he4.hydro.hst`, `he4.log`, `rt_profile.bin` and the driver's
own area-weighted `L_rad,out/L` face-budget print.  The `.bin`/`.rst`/`.cbin` were deleted
after measuring and the two job logs truncated to the startup block, the dt-collapse
blocks and the cycle diagnostics; every number here was read off the full log first.

Column key: `lnKE1` / `lnKEh` are ln of the hst radial and (2-KE + 3-KE) horizontal
kinetic energies; `dmass` is the relative drift of hst column 3; `floors` is
dfloor+efloor+tfloor+vceil cumulated from the event log (put on the time axis through the
driver's own `cycle=/time=` prints); `L/Lstar` is the pgen's exact area-weighted
`L_rad,out/L`; `d_rho_k` / `d_T_k` are the area-weighted shell means at the FeCZ base
(0.635 R = 1.5014e11, i = 15), the kappa peak (0.76 R = 1.8037e11, i = 29) and the FeCZ
top (0.97 R = 2.2985e11, i = 79), relative to t = 0; `T` is the pgen's `wtemp` slot, so
read the *relative* change only.  `vr_rms/vMLT` is sqrt(<v1^2>) at the kappa peak over
v_MLT = 1.45e7 cm/s.

## ARM A -- `inner_bc = open` -- FAIL at 0.117 turnover

```
  t/turn      t[s]     dt[s]    lnKE1   lnKEh   dmass    eos_fail floors  fofc   L/Lstar   d_rho_0 d_T_0  d_rho_1 d_T_1  d_rho_2 d_T_2   vr_rms/vMLT
   0.00         0      12.9   75.50   74.64  0.0e+00        0       0      0     0.938   0.000   0.000   0.000   0.000   0.000   0.000     0.000
   0.02     94.48      14.1   86.39   74.64  4.1e-04        0       0      0     0.786  -0.018   0.018  -0.000  -0.014  -0.001  -0.038     0.006
   0.04     195.5        14   88.14   74.64  1.1e-03        0       0      0     0.712   0.030   0.043  -0.017  -0.017  -0.004  -0.056     0.096
   0.06     289.8      10.5   90.45   74.90  2.1e-03        0       6      3     0.682   0.078   0.072  -0.032  -0.006  -0.010  -0.063     0.285
   0.08     376.8      8.62   91.64   82.71  3.5e-03        0      19     13     0.691   0.146   0.106  -0.021   0.010  -0.018  -0.062     0.296
   0.10     471.8      6.55   92.71   83.60  1.2e-02        0      24     13     0.728   0.008   0.149   0.041   0.038  -0.031  -0.056     0.058
   0.12       564   0.00551   95.59   95.35  2.1e-02        0   37763   4012   141.478  -0.128   0.152   0.138   0.079  -0.045   0.038     0.477
event-log totals: eos_fail=0 dfloor=99 efloor=37664 tfloor=0 vceil=0 fofc=4012 max c2p it=0
```

**Alive:** no.  `dt` falls off the 0.5 s gate between cycle 62 and cycle 87 and never
comes back (cycle 1000 at t = 566.4 s with dt = 2.5e-3 s).

**The first collapsing cell** (the driver's own dt-collapse print, rank 1):

```
### dt COLLAPSE cycle=62 time=548.852 dtold=4.80609 dt=1.12569 | hydro=1.12569
    hydro dt is set by cell (m,k,j,i) = (1,5,3,5) gid = 13
    r=1.22574e+11 rho=7.90847e-09 T=3.68382e+13 p=9.7397e+07 cs=1.28161e+08
    v=(3.31473e+08, -72497.3, -188870)
```

i = 5 is **six cells above `x1min` = 1.18585e11**, i.e. inside the 24-cell bottom sponge
and right on the open inner face.  The temperature is the garbage of a failed EOS
inversion, the radial velocity is +3.3e8 cm/s = 23 v_MLT outward and the transverse
components are 1e3 times smaller: **the runaway is purely radial and it starts at the
inner boundary.**  The next two collapses walk outward along the same column
(i = 33 at t = 560.6 s with v1 = 4.9e9, i = 50 at t = 561.9 s with v1 = 8.5e9).

**The term that runs away** is the open inner boundary's own mass and energy transfer.
`dmass` is **+2.1 % in 0.12 turnover** -- mass pours in through the inner face, while the
face budget shows energy *leaving* it (`inner gain/L` = -0.6e-2 .. -1.7e-2 throughout) and
the deep radiative luminosity at the two-stream cut climbing to `L_rad,cut/L` = 6.6.  The
radial kinetic energy grows **exponentially and linearly in ln**: lnKE1 goes 75.50 -> 92.71
in 472 s, an e-folding time of **27 s** = 0.006 turnover, i.e. ~170 e-folds per turnover.
Nothing convective grows at that rate; this is a boundary instability.  `eos_fail` stays 0
and the floors and FOFC are quiet (24 floors, 13 FOFC firings in 472 s) until the death
interval, in which 37 763 floor events and 4012 FOFC firings land at once -- FOFC is a
*symptom* of the collapse, not a precursor, and at 4012 firings over ~5.9e5 cells it is
never a spatially organised effect, so the cube-vertex question does not arise.  (No
`hydro_fofc` dump exists: the only `bin` written was t = 0.  If the spatial pattern ever
matters, run with `output4/variable=hydro_fofc output4/file_type=bin`.)

## ARM B -- `inner_bc = wall` + `rt_bottom_flux` -- FAIL at 0.366 turnover

```
  t/turn      t[s]     dt[s]    lnKE1   lnKEh   dmass    eos_fail floors  fofc   L/Lstar   d_rho_0 d_T_0  d_rho_1 d_T_1  d_rho_2 d_T_2   vr_rms/vMLT
   0.00         0      12.9   75.50   74.64  0.0e+00        0       0      0     0.938   0.000   0.000   0.000   0.000   0.000   0.000     0.000
   0.02     94.48      14.1   86.10   74.64  4.8e-08        0       0      0     0.786  -0.019   0.016  -0.000  -0.014  -0.001  -0.038     0.006
   0.04     196.3      14.8   87.86   74.64  1.3e-07        0       0      0     0.711  -0.060   0.028  -0.017  -0.022  -0.004  -0.057     0.073
   0.06     294.3        13   88.50   76.21  1.6e-07        0       0      0     0.673   0.126   0.036  -0.046  -0.022  -0.010  -0.066     0.246
   0.08     382.8      12.3   88.78   77.78  1.6e-07        0       0      0     0.662  -0.034   0.042  -0.070  -0.021  -0.019  -0.070     0.361
   0.10     478.3        12   88.87   79.65  9.0e-08        0       0      0     0.660   0.219   0.047  -0.097  -0.020  -0.032  -0.073     0.387
   0.12     570.4      11.3   89.03   80.16 -2.7e-08        0       0      0     0.665   0.060   0.051  -0.130  -0.019  -0.048  -0.076     0.311
   0.14     659.6      11.1   88.92   80.32 -1.5e-07        0       1      0     0.673  -0.041   0.053  -0.154  -0.017  -0.066  -0.079     0.162
   0.16     758.4      10.9   88.62   80.52 -2.8e-07        0       3      0     0.680   0.148   0.057  -0.158  -0.015  -0.089  -0.084     0.006
   0.18     854.3      10.5   88.97   80.90 -4.1e-07        0       7      0     0.686   0.603   0.059  -0.135  -0.014  -0.114  -0.091     0.076
   0.20     949.6      10.3   88.80   81.75 -5.7e-07        0       8      0     0.689   0.495   0.060  -0.132  -0.012  -0.142  -0.098     0.079
   0.22      1040      9.94   88.90   81.72 -7.3e-07        0      10      0     0.694   0.371   0.061  -0.153  -0.010  -0.172  -0.106     0.079
   0.24      1130      9.78   89.00   82.18 -8.8e-07        0      14      0     0.699   0.264   0.061  -0.209  -0.009  -0.204  -0.115     0.075
   0.26      1227      9.56   89.07   82.68 -1.0e-06        0      33      0     0.702   0.117   0.063  -0.276  -0.007  -0.242  -0.125     0.154
   0.28      1323      8.28   89.49   82.94 -1.1e-06        0     111      0     0.704   0.135   0.066  -0.275  -0.005  -0.280  -0.136     0.229
   0.30      1413      6.96   89.58   83.09 -1.2e-06        0     286      0     0.706   0.271   0.068  -0.200  -0.003  -0.318  -0.147     0.223
   0.32      1506       5.8   89.60   83.23 -1.3e-06        0     740      0     0.712   0.384   0.069  -0.058   0.001  -0.357  -0.158     0.076
   0.34      1603      5.91   89.65   83.40 -1.4e-06        0    1089      0     0.718   0.483   0.068   0.159   0.003  -0.399  -0.170     0.349
   0.36      1693         4   89.74   83.59 -1.5e-06        0    1230      0     0.721   0.280   0.066   0.150   0.004  -0.436  -0.182     0.349
event-log totals: eos_fail=0 dfloor=0 efloor=1230 tfloor=0 vceil=0 fofc=0 max c2p it=0
```

**Alive:** to t = 1723 s (0.366 turnover), then no.  `dt` holds 10-15 s for a third of a
turnover -- comfortably above the 0.5 s gate -- decays to 4 s over the last 0.1 turnover
and collapses to 1.3e-1 s at cycle 190, 2e-2 at cycle 205 and 4.6e-4 by cycle 2215.

**The first collapsing cell:**

```
### dt COLLAPSE cycle=190 time=1723.03 dtold=1.26915 dt=0.130482 | hydro=0.130482
    hydro dt is set by cell (m,k,j,i) = (1,13,15,24) gid = 13
    r=1.63758e+11 rho=5.37469e-10 T=5.24607e+11 p=69.4262 cs=465713
    v=(-5.08928e+09, 5.28677e+06, -2.3981e+06)
```

then, three seconds later,

```
### dt COLLAPSE cycle=205 time=1726.11 ... cell (1,11,16,54) gid = 13
    r=2.1102e+11 rho=9.72695e-09 T=1.56485e+14 p=3.16444e+10 cs=2.08272e+09
    v=(1.17021e+10, 1009.81, -430.196)
```

**This is the same death as the 1-D gate**, at the same place.  `../README.md` section 2b
records `(1,3,4,53)`, `r = 2.10e11`, purely radial `v1 = 8.8e9` and the same
`T = 2.63e14` EOS-inversion sentinel; here it is `i = 54`, `r = 2.11e11`, `v1 = 1.17e10`
with transverse components 1e7 times smaller, and the sentinel reappears verbatim at the
two later collapses.  Going to 3-D bought **0.37 turnover instead of 0.47** -- i.e.
nothing.

**The term that runs away** is the **energy floor in an emptying upper envelope**, and the
run-up is visible for a fifth of a turnover before the death:

* `efloor` is the only counter that moves.  It is 0 through t = 570 s, then
  1, 3, 7, 8, 10, 14, 33, 111, 286, 740, 1089, 1230 -- doubling roughly every 0.02
  turnover.  `eos_fail`, `dfloor`, `tfloor`, `vceil` and `fofc` are **all exactly zero**
  for the whole run, and `max c2p it` never leaves 0.
* The shell means say why.  The **FeCZ top drains monotonically from t = 0**: at 0.97 R,
  rho **-43.6 %** and T **-18.2 %** by 0.36 turnover (e-folding ~0.6 turnover), with
  `<v1> = -2.6e6` cm/s everywhere -- the outer envelope is falling inward at 0.18 v_MLT.
  The FeCZ base and kappa peak take it up: rho +28 % and +15 %.  The mass drift is
  -1.5e-6, so nothing leaves the domain; the column is simply restratifying downward.
* The emergent luminosity **never recovers**: `L_rad,out/L` falls 0.938 -> **0.660** in the
  first 0.10 turnover and limps back only to 0.721.  A third of `L` is being absorbed and
  is not re-radiated, which is exactly the 1-D diagnosis -- the FeCZ is heated from below
  faster than a non-convecting column can radiate.
* The floored cells then carry a garbage inverted temperature, their pressure collapses
  (p = 69 at rho = 5.4e-10), and the neighbours pour in at 5e9 cm/s.

**Convection never starts.**  `lnKEh` rises only 74.6 -> 83.6, i.e. KE_horizontal reaches
3.1e36 against KE_radial 9.4e38 -- **0.3 %**.  `vr_rms` at the kappa peak oscillates
between 0.08 and 0.39 v_MLT with a period of ~0.1 turnover: that is a coherent radial
pulsation of the whole column (the same signature the 1-D run had), not overturning.  A
1e-3 c_s seed needs several turnovers to reach v_MLT; the run has 0.37.

## Verdict

**Use the WALL inner boundary with `rt_bottom_flux = true` and
`rad_flux_inner = 1.305278e15` (arm B).**  It is better on every measure:

| | A (open) | B (wall + base flux) |
| --- | --- | --- |
| lifetime | 0.117 turnover | **0.366 turnover** (3.1x) |
| mass drift to the death | **+2.1 %** | -1.5e-6 |
| ln KE1 e-folding | 27 s (170 /turnover) | no exponential; KE1 flat from 0.04 turnover |
| first collapsing cell | i = 5, ON the inner face | i = 24 / i = 54, in the FeCZ |
| floors before the death | quiet, then 3.8e4 at once | efloor doubling for 0.2 turnover |
| failure mode | its own inner-boundary runaway | the 1-D failure, unchanged |

Arm A does not merely fail; it adds a violent, purely radial, inner-boundary-driven
instability **on top of** the failure arm B shows, and it loses mass conservation while
doing it.  As configured (`s_relax_cs = 0.1`, `s_relax_cp = 0.3`, `rad_flux_inner = 0`,
24-cell bottom sponge) the deep-adiabat relaxation is not usable for this star: it is a
net energy *sink* (`inner gain/L` < 0 at every sample) while being a mass *source*.

Arm B's inner boundary is therefore innocent of the death, and the death is the one the
1-D gate already identified.  **The inner boundary is not what is blocking this run.**

## Recommended next step

**Another smoke iteration, not the production grid.**  The production grid (nx2 = nx3 =
320) would cost ~100x and reproduce the same collapse at the same 0.37 turnover; nothing
in these tables is resolution-limited -- `eos_fail = 0`, `fofc = 0`, `max c2p it = 0`, the
cube vertices are silent, and the failure is a one-dimensional thermal one.

Two things must change, in this order:

1. **Rebuild the initial column so the upper envelope can hold itself up.**  The
   diagnostic number is `L_rad,out/L` = 0.66 by 0.06 turnover with rho at 0.97 R falling
   44 % in a third of a turnover, *before* any convective flux exists.  This is
   `../README.md` step-2 option 1: relax against a radiative-equilibrium target (push
   `ic_tau_rad` through the domain, or relax the 1-D column with the wall + base flux and
   hand the result on with `../mk_ic_from_profile.py`, which was written for exactly this
   and is still unused).  The 3-D run then starts from a state it can hold while the
   convection grows into it, instead of racing a 0.6-turnover drain.
2. **Only then raise the seed** to `vpert = 1e-2` (the arm-C variant).  A stronger seed
   alone cannot help: at 0.37 turnover of lifetime even v = v_MLT at t = 0 would not have
   time to build a convective flux, and the drain starts at t = 0 regardless of the seed.

A cheap intermediate worth one 15-minute `apudev` job once (1) is in hand: arm B again on
this same smoke grid, gating only on whether `L_rad,out/L` stays within 5 % of 1 and
`efloor` stays at 0 through one turnover.  If it does, run arm C (vpert 1e-2) for three
turnovers and look for `lnKEh` catching `lnKE1`; only then go to the production grid.

## Files

* `chain.sh` -- the submit script (restart-chaining, `ARM`/`WALL` through `--export`).
* `gate.py <armdir> [label]` -- the table generator.
* `A/`, `B/` -- `he4.hydro.hst`, `he4.log` (the event log), `rt_profile.bin`,
  `rt_surface.bin`, `column_he4_<arm>.txt`.  `.bin`/`.rst`/`.cbin` deleted.
* `tableA.txt`, `tableB.txt` -- the tables above as generated.
* `A.11758720.log`, `B.11758721.log` -- the job logs, truncated (startup block, all
  dt-collapse blocks, the cycle diagnostics).
