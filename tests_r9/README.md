# he4_presn: THE 0.97 R PILE-UP IS THE GREY TWO-STREAM'S OWN FLUX ERROR, NOT THE TAPER,
# NOT THE CLOSURE AND NOT A BOUNDARY

2026-09-17/18, viper, branch `he4-presn-global`, against `a6d66c3c` (= the `tests_r8`
commit).  Read `tests_r8/README.md` sections 4.1, 6 and 8 first: this round answers its
section 7 item 4 (`mlt_alpha = 0` / `rt_bottom_flux = false` and the four other bisection
arms) and its section 4.2 (the `rt_profile.bin` `eint` slot), and collects `x_gate5`.

Binaries: `tests_r9/athena_v5` (HIP, `Kokkos_ARCH_AMD_GFX942_APU`, MPI,
`PROBLEM=red_giant`; slot 6 fixed) and `tests_r9/athena_v6` (the same with slot 5 fixed
too -- the one to keep).  Every arm is 1-D: `nx2 = nx3 = 4`, `meshblock 96x4x4`,
`vpert = 0`, the `tests_r8` `x_gate` overrides (`inner_bc = wall`,
`rt_bottom_flux = true`, `rad_flux_inner = 1.305278e15`, `mlt_alpha = 1.5`,
`rt_rad_force = true`, the gated taper), one rank, one MI300A.  Turnover = 4705 s.
Every `.bin`, `.rst` and `.cbin` has been deleted after measuring; what is kept is the
`.hst`, the event log, `column_*.txt`, `mltfaces_*.txt`, `rt_profile.bin` and
`rt_surface.bin` (a few hundred kB per arm).

---

## 0. THE ANSWER IN FOUR NUMBERS

| | measured |
| --- | --- |
| the star's own EMERGENT luminosity, once the two-stream has relaxed (0.1 turnover) | **L_out = 0.981 L**, falling to **0.956 L** by 1.4 turnovers |
| the thermal time of the whole domain, E_domain/L | **0.462 turnover** (E_0 = 5.016e41 erg, L x turnover = 1.085e42 erg) |
| the RADIATIVE COOLING TIME of the 0.99 R cell at t = 0, from the run's own face dump | **6.5 s**, against `dt` = **12.9 s** at cycle 0 |
| `F_2s/F_req` at t = 0, from 0.98 R to the top | **1.007 -> 1.042 -> 1.21 -> 1.28** |

A grey atmosphere must carry `F = L/(4 pi r^2)` at EVERY radius.  This column does, to
0.1 %, from 0.51 R to 0.94 R.  Then the two-stream carries 28 % too much at the top face
at t = 0 and, after it has relaxed against the gas, 2-4 %
too little out of the star.  The domain's thermal time is half a turnover, so a 2-4 %
shortfall is a few per cent of the envelope's own energy per turnover: the upper envelope
heats, the iron opacity rises with it, the local Eddington ratio at 0.94 R crosses 1.0 at
one turnover, and the layer lifts.  That is the "0.97 R pile-up".  **It is numerical, the
operator is the grey spherical two-stream, and NONE of the six bisection arms removes
it** -- because none of them changes that operator's flux.

---

## 1. `x_gate5`: THE 5-TURNOVER ARM DID NOT SETTLE, AND IT DID NOT DIE EITHER

`tests_r8/r8_1d_long.sh`, arm `x_gate5/`, re-submitted 22:53 for a 3 h `apu` slot.

| | |
| --- | --- |
| reached | **t = 11844 s = 2.518 turnovers** (of the 5 asked for) |
| `dt` at the end | **0.0757 s**, from 12.93 at cycle 0 -- a factor 171 |
| `dt COLLAPSE` print | **none** |
| `eos_fail` | **0** |
| domain mass | 1.9125e26 -> 1.8824e26 g = **-1.57 %** |
| total energy | 5.0158e41 -> 6.019e41 erg = **+20 %** |
| `eos_dfloor` / `eos_efloor` per 0.01 turnover, at the end | 3.85e6 / 1.19e6 |

So the 5-turnover settling gate **FAILS**, for the reason `tests_r8` section 8 predicted:
not wall time, the 0.97 R layer.  At 0.0757 s a turnover is 6.2e4 cycles.  The energy
budget is the tell: **+20 % of the domain's energy in 2.5 turnovers**, which at
L_out/L = 0.96-0.98 is exactly the accumulated flux shortfall (below).

## 2. WHAT IS PILING UP, WHERE, AND WHY -- THE CHARACTERISATION

All of section 2 is the `x_gate` configuration: `x_gate5` for the history (53 profile
records over 2.6 turnovers) and `base9` for the t = 0 face budget (identical run, new
binary; the two agree to 1.7e-4 in mass and 0.18 in `dt` by 2.2 turnovers, i.e. they are
the same run diverging in the chaotic phase after 1.5 turnovers).

### 2.1 THE t = 0 FACE BUDGET, FROM THE RUN'S OWN `mltfaces` DUMP

`F_cond` is **exactly zero at every face** -- `w_blend = 0.000` over the whole column,
which is what `rad_tau_lo/hi = 1e5/1e6` (the shipped whole-column blend) means: the
radial conduction operator is inert and the two-stream owns all 96 cells.  The MLT
closure supplies 10.1 % of L at 0.90 R and **nothing at all above 0.94 R**, because the
two-stream is already over-carrying there and the deficit `max(0, F_req - F_cond - F_res
- (1-w)F_2s)` is zero.  `F_res` is 1e-11 of L (it is a 1-D column).

| face r/R | `F_req` [cgs] | `F_2s/F_req` | `F_cond/F_req` | `F_MLT/F_req` | `L_face/L` |
| --- | --- | --- | --- | --- | --- |
| 0.899 | 4.038e14 | 0.899 | 0.000 | 0.101 | 1.0000 |
| 0.940 | 3.695e14 | 1.001 | 0.000 | 0.000 | 1.0011 |
| 0.960 | 3.541e14 | 1.003 | 0.000 | 0.000 | 1.0031 |
| 0.970 | 3.465e14 | 1.003 | 0.000 | 0.000 | 1.0029 |
| 0.981 | 3.388e14 | 1.007 | 0.000 | 0.000 | 1.0066 |
| 0.990 | 3.331e14 | 1.042 | 0.000 | 0.000 | **1.0422** |
| 1.001 | 3.257e14 | 1.213 | 0.000 | 0.000 | **1.2125** |
| 1.009 | 3.204e14 | 1.278 | 0.000 | 0.000 | **1.2781** |

And DEEPER, which is the half of the budget that is in fine shape and has to be said so:

| face r/R | 0.513 | 0.599 | 0.638 | 0.704 | 0.757 | 0.802 | 0.851 | 0.899 | 0.921 | 0.940 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `F_2s/F_req` | 1.000 | 1.000 | 0.998 | 0.921 | 0.873 | 0.852 | 0.848 | 0.899 | 0.976 | 1.001 |
| `F_MLT/F_req` | 0.000 | 0.000 | 0.002 | 0.079 | 0.127 | 0.148 | 0.152 | 0.101 | 0.024 | 0.000 |
| **SUM** | **1.000** | **1.000** | **1.000** | **1.000** | **1.000** | **1.000** | **1.000** | **1.000** | **1.000** | **1.001** |

**The t = 0 budget is EXACT from 0.51 R to 0.94 R** -- the whole-column two-stream carries
0.85-1.00 of L and the shell-mean MLT closure supplies the 0-15 % residual in exactly the
Schwarzschild-unstable shells (`grad - grad_ad` turns positive at 0.638 R, and that is
where `F_MLT` turns on).  The `tests_3d/handover` complaint that the whole-column blend
"carried only 0.10-0.74 of L over 0.55-0.90 R" does NOT apply to this column: the inner
wall was moved out to 0.50 R since, and the hole is gone.  **The entire t = 0 defect is
in the last five cells**, above tau = 3.

### 2.2 THE DIVERGENCE OF THAT, PER SHELL: THE TOP TWO CELLS COOL FASTER THAN `dt`

Between the cell's OWN two faces, at t = 0, with the cell's own thermal content
`eint x V`:

| cell r/R | its faces r/R | V [cm^3] | `eint` [erg/cm^3] | E_shell [erg] | L_in - L_out [erg/s] | /L | cooling(-)/heating(+) time |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.9410 | 0.9397-0.9423 | 3.77e32 | 1.761e6 | 6.63e38 | **+8.98e34** | +0.0004 | heats in **1.57 turnover** |
| 0.9613 | 0.9600-0.9626 | 4.00e32 | 7.834e5 | 3.14e38 | **+4.82e34** | +0.0002 | heats in **1.38 turnover** |
| 0.9718 | 0.9705-0.9732 | 4.25e32 | 4.278e5 | 1.82e38 | **+4.30e34** | +0.0002 | heats in **0.90 turnover** |
| 0.9827 | 0.9814-0.9841 | 4.49e32 | 2.036e5 | 9.15e37 | **-1.204e36** | -0.0052 | cools in **0.016 turnover = 76 s** |
| 0.9911 | 0.9897-0.9925 | 4.63e32 | 9.493e4 | 4.40e37 | **-6.761e36** | -0.0293 | cools in **0.0013 turnover = 6.5 s** |
| 1.0023 | 1.0010-1.0037 | 4.66e32 | 1.510e4 | 7.04e36 | **-1.028e37** | -0.0044 | cools in **0.0007 turnover = 3.4 s** |

`dt` at cycle 0 is **12.93 s**.  The cells at 0.99 R and 1.00 R are handed a radiative
sink that would empty them in HALF a timestep and a QUARTER of one.  They are held only
by the implicit column solve (`rt_implicit_column = 3`), which is why nothing NaNs -- but
they cannot stay in hydrostatic balance, and they do not.

### 2.3 THE HISTORY: THE 0.99 R CELL EMPTIES FIRST, THEN THE WHOLE ENVELOPE LIFTS

`rho/rho(t=0)` per 0.25 turnover, `x_gate5` (`tests_r9/prof9.py`, which corrects both
pre-fix `rt_profile` slots, see section 4):

| t/turn | 0.941 R | 0.959 R | 0.969 R | 0.980 R | 0.991 R | 1.000 R | 1.010 R |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.25 | 1.05 | 0.99 | 1.01 | 0.97 | 1.04 | 1.13 | **0.58** |
| 0.50 | 1.01 | 1.01 | 0.97 | 0.96 | 1.03 | 1.09 | **0.41** |
| 0.75 | 1.06 | 1.06 | 1.11 | 0.91 | 0.94 | 1.05 | **0.40** |
| 1.00 | 1.22 | 1.04 | 1.07 | 1.27 | **0.37** | 1.39 | 0.52 |
| 1.25 | 1.06 | 1.66 | 1.34 | 1.03 | 1.08 | 1.62 | 0.57 |
| 1.50 | **0.68** | 1.44 | **2.36** | **3.32** | 0.92 | 1.79 | 1.92 |
| 1.75 | **0.42** | 0.84 | 1.69 | **3.39** | **6.09** | **8.60** | **20.8** |
| 2.00 | 0.46 | 0.88 | 1.48 | 3.06 | 5.94 | 7.78 | 10.3 |
| 2.25 | 0.31 | 0.18 | 0.18 | 0.24 | 0.32 | 1.85 | 0.001 |

and the same cells' `T`, `v1`, `Gamma = kappa F_2s/(c g)` and taper `w`:

| t/turn | T(0.941R) | Gamma(0.941R) | v1(0.980R) [cm/s] | Mdot(0.980R) [g/s] | eint(0.969R) | w(1.010R) |
| --- | --- | --- | --- | --- | --- | --- |
| 0.00 | 1.188e5 | 0.930 | 0 | 0 | 5.04e5 | 0.000 |
| 0.50 | 1.195e5 | 0.948 | -1.02e5 | -1.9e20 | 4.84e5 | 0.000 |
| 1.00 | 1.238e5 | **1.039** | +7.10e5 | +1.7e21 | 5.42e5 | 0.000 |
| 1.25 | 1.323e5 | **1.095** | -5.9e3 | -1.3e19 | 6.70e5 | 0.000 |
| 1.50 | 1.389e5 | 1.013 | **+3.06e6** | **+1.96e22** | 1.43e6 | 0.000 |
| 1.75 | 1.403e5 | 0.821 | **+3.24e6** | **+2.12e22** | 2.13e6 | **1.000** |
| 2.00 | **1.983e5** | 0.498 | **-7.55e6** | **-4.46e22** | **1.27e7** | 1.000 |
| 2.25 | 1.160e5 | 4.15 | -4.82e6 | -2.2e21 | 4.62e5 | 1.000 |

**PRECISELY WHAT HAPPENS, IN ORDER.**

1. **0 to 0.5 turnover -- the top cell drains.**  The 1.010 R cell, whose t = 0 radiative
   cooling time is 3.4 s, loses 60 % of its density in a quarter turnover.  This is the
   28 % over-carry of section 2.1 seen from the inside; from the outside it is
   `L_out/L = 1.278` at the t = 0 dump collapsing to **0.981 at the 0.1-turnover dump**
   (section 3).  The mass it loses goes to the cells below.
2. **0.5 to 1.25 turnovers -- the 0.99 R cell empties, and the interior heats.**  The
   0.991 R cell -- the one with the 6.5 s cooling time -- is at 0.37 of its initial
   density at one turnover, and the 0.980 R cell below it is at 1.27.  That is the
   pile-up, and it is **between the faces at r/R = 0.9814 and 0.9841 (2.3271e11 and
   2.3336e11 cm)**: the cell that gains is the one immediately BELOW the cell with the
   largest flux divergence.  Meanwhile `T` at 0.94 R rises 1.188e5 -> 1.238e5 -> 1.323e5,
   the Rosseland opacity rises with it (`tau` at that cell goes 42 -> 50 -> 66), and
   **`Gamma` at 0.941 R crosses 1.0 at exactly one turnover** (0.930 -> 1.039 -> 1.095).
3. **1.25 to 1.75 turnovers -- the envelope LIFTS.**  Every cell from 0.94 to 1.01 R
   acquires `v1 = +2` to `+4.3e6 cm/s` (outward; `v_MLT` = 1.45e7, `c_s` ~ 3e6 there, so
   this is sonic), `Mdot(0.98 R) = +2.1e22 g/s`, `rho` reaches 3-21x its initial value
   above 0.98 R, `eint` at 0.969 R reaches 4.2x, and `T` inverts (T rises OUTWARD).  The
   pile MOVES OUT: its peak is at 0.980 R at 1.5 turnovers and above 1.01 R at 1.75.
4. **2.0 turnovers -- it falls back.**  The whole column goes isothermal at 2.0e5 K with
   `eint` 25-140x its initial value, `v1 = -7.5e6 cm/s` at every radius,
   `Mdot = -4.5e22 g/s`, `F_2s` REVERSES SIGN (the T inversion), and `dt` falls to 0.1 s.
   This is the event `tests_r8` section 4.1 found between 1.75 and 2.0 turnovers.

**Cooling (contracting) or pushed (force)?  PUSHED, and HEATED.**  Over 1.25-1.75
turnovers the piling shells have `v1 > 0` (outward), `Mdot > 0`, `eint` rising 2-30x and
`Gamma` at the layer below them at or above 1.0.  Nothing in the pile-up phase is a
cooling contraction.  The one cell that IS cooled -- 0.99-1.01 R, the 3-7 s sink -- drains
and its mass lands below; that is the seed, not the pile.

### 2.4 THE ENERGY BUDGET OF THE PILE-UP SHELL, PER TURNOVER

For the 0.9718 R cell (V = 4.25e32 cm^3, E_shell(0) = 1.82e38 erg), in erg/cm^3/s:

| term | value | how |
| --- | --- | --- |
| radiative deposit at t = 0, `-div F_carried` | **+101** | the run's own `mltfaces` faces, section 2.2 |
| MLT closure divergence | **0** | `F_MLT = 0` at both faces (the deficit is zero there) |
| conduction divergence | **0** | `w_blend = 0`, the operator is inert |
| resolved enthalpy-flux divergence | < 1e-6 of the above | 1-D, `F_res/F_req = 1e-11` |
| **measured** `d eint/dt`, 0 -> 0.25 turnover | **-4.3** | `rt_profile` slot 6 |
| **measured** `d eint/dt`, 1.25 -> 1.50 turnover | **+651** | ditto |
| **measured** `d eint/dt`, 1.75 -> 2.00 turnover | **+9025** | ditto |

So at t = 0 the shell's only energy source is the +101 erg/cm^3/s radiative deposit --
0.02 % of L, but 0.90 turnover of its own thermal content -- and the early measured drift
(-4.3) is the PdV of the settling, an order below it.  By 1.5 turnovers the measured
heating is 6x the t = 0 radiative deposit and by 2.0 turnovers 90x: past one turnover the
shell is no longer being heated by its own flux divergence but by the compression of the
falling envelope.  **The radiative deposit is the trigger, the dynamics are the event.**

## 3. THE EMERGENT LUMINOSITY, WHICH IS THE WHOLE STORY

`tests_r9/lout9.py` on `rt_surface.bin` (area weighted over the 96 angular columns; in
1-D they agree to 1e-3 with each other until 1.5 turnovers):

| t/turn | `base9` | `a_mlt0` | `c1_topwall` | `d_taperup` |
| --- | --- | --- | --- | --- |
| 0.00 (the ic, before the solver has run) | **1.2779** | 1.2779 | 1.2779 | 1.2779 |
| 0.10 | **0.9813** | 0.981 | 0.985 | **0.9945** |
| 0.50 | 0.9775 | 0.978 | 0.978 | 0.9917 |
| 0.90 | 0.9702 | 0.968 | 0.978 | 0.9723 |
| 1.30 | 0.9599 | 0.930 | 0.938 | 0.9616 |
| 1.40 | 0.9563 | -- | -- | -- |
| 1.70 | (1.5: 0.9696) | 0.0010 | 1.019 | 1.027 |

Two separate defects are visible here.

* **The 28 % startup jump.**  The initial column is a grey EDDINGTON atmosphere
  (`f = 1/3`, shot so that `tau = 2/3` sits at R -- `column_he4_presn_sph.txt` header),
  and the solver is a DISCRETE-ORDINATE two-stream with `ck_nquad = 2`
  (`mu = +-1/sqrt(3)`).  The two closures do not have the same `T(tau)` relation, so the
  solver's emergent flux on the ic state is 1.278 L, and it relaxes to its own solution
  within one dump interval -- kicking the 3.4 s and 6.5 s cells on the way.  This is an
  **ic/solver inconsistency**, and it is exactly the `L_out/L` step the arms all share.
* **The 2-4 % standing shortfall.**  After relaxation the star radiates 0.98 L and decays
  to 0.956 L.  This is the documented spherical-form error of the grey two-stream:
  `he4_presn_cs.athinput` states, for the whole-column blend, "above tau = 20, where the
  two-stream owns the column, `1 - H_T/(2r)` is 0.96-0.97, so the spherical-form error
  there is 3-4 %".  **The measurement is 0.956-0.981.**  With the domain thermal time at
  0.462 turnover, 2-4 % of L is 4-9 % of the envelope's own energy per turnover, and the
  `x_gate5` history confirms it: **+20 % of the domain energy by 2.5 turnovers**.

### 3.1 AND IT IS NOT THE TAPER -- THE ONE CLEAN DISCRIMINATION

In the baseline the over-carry looks like a taper effect: `F_2s/F_req - 1` fits
**0.262 (1-w) + 0.004 with correlation 0.989** over the nine faces where `w < 1`, and
saturates at +0.278 where `w = 0`.  But the taper window was READ OFF tau 3 and 0.3, so
`w` and `tau` are degenerate in that fit.  Arm (d) breaks the degeneracy by moving the
window to tau 1 / 0.1, i.e. 0.997-1.009 R instead of 0.986-1.005 R:

| face r/R | `w` (baseline) | `w` (arm d) | `F_2s/F_req - 1` (baseline) | `F_2s/F_req - 1` (arm d) |
| --- | --- | --- | --- | --- |
| 0.9869 | 0.985 | **1.000** | +0.02280 | **+0.02279** |
| 0.9925 | 0.754 | **1.000** | +0.07147 | **+0.07147** |
| 0.9982 | 0.328 | **0.964** | +0.15888 | **+0.15890** |
| 1.0010 | 0.135 | **0.821** | +0.21252 | **+0.21254** |
| 1.0065 | 0.000 | **0.221** | +0.27610 | **+0.27613** |
| 1.0118 | 0.000 | 0.000 | +0.27797 | **+0.27801** |

**The flux profile does not move at all** -- five significant figures, at faces where `w`
changed from 0.13 to 0.82.  The over-carry is a function of the RADIUS (of `tau`), not of
the taper weight.  The `0.262 (1-w)` fit was a coincidence of how the window was chosen.
The taper, the EOS gate and `rt_rad_force` are all exonerated for this defect.

## 4. TASK 3: THE `rt_profile` SLOTS.  BOTH 5 AND 6 WERE WRONG

`tests_r8` section 4.2 left open that slot 6 is "1.02x the ic file's `eint` at 0.50 R and
80x at 1.01 R".  Both it and slot 5 are now fixed in `red_giant`'s own writer
(`src/pgen/red_giant.cpp`, `RedGiantProfileDump`); the box's writer is untouched.

| slot | what it WAS | what it IS | why it was wrong |
| --- | --- | --- | --- |
| 5 | `wtemp` RAW = the CODE temperature | `temp_cgs * wtemp` = **kelvin** | `Hydro::wtemp` is `k_B T/(mu_ref m_u)` in code energy/mass; `EOS_Data::temp_cgs` (`eos.cpp:281`, `= temperature_cgs()/mu`) is the kelvin per unit of it, `m_u/k_B = 1.2027e-8` on a cgs `<units>` block.  Raw, slot 5 was **8.3145e7 times** the kelvin its comment promised: T = 1.19e5 K at 0.94 R read out as 9.88e12. |
| 6 | `u0(IEN) - KE` | `EintFromCons(...)` = `u0(IEN) - KE(metric) - emag - rho Phi` | `<hydro>/etotgrav = true` in every `he4_presn` input, so `u0(IEN)` carries `rho Phi` with `Phi = GM(1/rin - 1/r) >= 0` (`red_giant.cpp:707`).  Slot 6 was `eint + rho Phi`: right AT the inner wall where `Phi = 0`, **80x** the true `eint` at the top of the domain, where `eint` is 1e-8 of the base and `rho Phi` is not. |
| 7 | `v1*(eint + w0(IPR))` | `v1*(eint + wder(IDPR))` | `IPR == IEN == 4`, and under a GENERAL EOS the primitive slot holds the INTERNAL ENERGY, not the pressure (`general_c2p_hyd.hpp:56`, `w.e = u.e - e_k`; the pressure is in `wder(IDPR)`).  Slot 7 was `v1*2*eint`, not the enthalpy flux. |

The record layout is documented in full at `red_giant.cpp:496-536`.

### 4.1 THE READ-BACK, ON THE FIXED BINARY

Arm `slotchk/` (`athena_v6`, 12 min on `apudev`), cycle-0 record against the ic file and
the run's own `column_slotchk.txt`:

| quantity | 0.50-0.95 R | worst in the domain |
| --- | --- | --- |
| slot 0 `rho` / ic `rho` | 0.99993-1.00000 | 0.99993 |
| slot 5 `T` / column `T` | **1.00017 .. 0.99855** | **0.9166** at r/R = 1.0051 |
| slot 6 `eint` / ic `eint` | **1.00056 .. 0.99518** | **0.8895** at r/R = 0.9968 |
| slot 7 / (`v1` x slot 6) at 0.1 turnover | **1.315** | -- |

Slot 7's 1.315 is the check that matters for that slot: `1 + p/eint = 1.315` means
`p/eint = 0.315`, i.e. 1/3, which is a radiation-dominated gas -- against the **2.000**
the old arithmetic gave by construction.  Slots 5 and 6 now reproduce the column and the
ic to 0.05-0.5 % over 0.50-0.95 R.  The residual 8 % in T and 11 % in `eint` is confined
to 0.99-1.01 R, the taper window, where the offline bilinear solve of
`make_ic_sph.py` and the run's bicubic-Hermite tapered EOS differ most -- the same window
section 3 is about, and one more reason not to trust anything in it to better than 10 %.

### 4.2 `an6.py` WAS NEVER WRONG, AND HERE IS WHY

`tests_r6/an6.py` reads slots 0, 1, 2, 3 and 5 and **never touches 6 or 7**, and it uses
slot 5 only through `fac = T_column/Q[0,5,:]`, i.e. it CALIBRATES the slot against the
initial column at t = 0 and reports `fac * Q[j,5,i]`.  That ratio is unit-free, so its
`T[K]` and `Gamma` columns were correct with the raw code temperature and are bitwise the
same kind of correct now (`fac` simply changes from 1.2e-8 to 1.0003).  Re-run on
`slotchk` and on `base9` it gives the same T and Gamma as `prof9.py`.  **The consumer that
was actually wrong is any reader that takes slot 5 or 6 at face value** -- `icchk.py`
(which is how the defect was found) and the new `prof9.py`, which is why `prof9.py` has a
`--raw` / `--traw` switch that corrects pre-fix files from either binary.

## 5. TASK 2: THE BISECTION.  SIX ARMS, AND NOT ONE OF THEM REMOVES THE PILE-UP

`r9_arm.sh`, one 1-D arm per 1-GPU `apu` job, `tlim` = 3 turnovers, a 1 h wall.  Every
number is at the SAME simulation time in every arm, which is what makes them comparable;
the arms that are still short of 3 turnovers were cut by the wall, not by physics (no arm
printed `dt COLLAPSE` and every arm has `eos_fail = 0`).

| arm | what it changes | reach [turn] | `dt` end | rho(0.97R)/rho0 at 1 / 1.5 / 2 turn | `dt` at 0.5 / 1 / 2 turn | `L_out/L` at 0.5 / 1 turn | M/M0 end | VERDICT on the pile-up |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `base9` | -- (the `x_gate` baseline, new binary) | 2.36 | 0.078 | **1.07 / 2.36 / 1.47** | 12.4 / 10.8 / 2.52 | 0.978 / 0.970 | 0.984 | -- |
| (a) `a_mlt0` | `mlt_alpha = 0` | **3.00** | 0.807 | **0.44 / 0.82 / 1.04** | 11.3 / 7.24 / 3.39 | 0.988 / **0.888** | **0.786** | **REPLACED by draining**: no pile, but the shell is at 0.44 at one turnover and 21 % of the star has left |
| (b) `b_botflux` | `rt_bottom_flux = false`, `rad_flux_inner` unchanged | 1.58 | 0.42 | **0.04 / 0.06 / --** | 15.5 / 2.39 / -- | **0.094** / 0.387 | 1.000 | **INVALID ARM, and a real finding**: see 5.1 |
| (c1) `c1_topwall` | `outer_bc = wall` | **3.00** | 2.84 | **1.07 / 2.38 / 1.01** | 12.4 / 10.8 / 4.64 | 0.978 / 0.970 | **0.664** | **LEFT EXACTLY WHERE IT WAS** (2.38 against 2.36) |
| (c2) `c2_nosponge` | `vdamp_top_tau = 0` | 1.60 | 1.25 | **0.96 / 1.85 / --** | **0.561** / 0.997 / -- | **1.003 / 0.993** | 0.997 | **LEFT** (1.85 against 2.36); best `L_out`, worst `dt` |
| (d) `d_taperup` | taper window moved to tau 1 / 0.1 (+ a rebuilt ic) | 2.35 | 0.083 | **1.09 / 2.38 / 1.42** | 12.4 / 10.7 / 2.85 | 0.992 / 0.973 | 0.978 | **LEFT, to five significant figures** -- see 3.1 |
| (e) `e_handover` | `rad_tau_lo/hi = 20/300` | 1.18 | 0.226 | **0.04 / -- / --** | 14.4 / **0.235** / -- | **0.024** / 0.970 | 1.000 | **REPLACED by draining**, and it dies first |
| (f) `f_nx144` | `nx1 = 144` (the same stretch map) | 2.26 | 0.061 | **0.99 / 2.39 / 1.63** | 8.33 / 7.34 / 1.51 | 0.981 / 0.960 | 0.987 | **LEFT, unchanged at 1.5x the radial resolution** |

`T(0.97R)/T0` at 2 turnovers is 2.40 (`base9`), 2.40 (`a_mlt0`), 2.36 (`d_taperup`), 2.43
(`f_nx144`) and 1.68 (`c1_topwall`): the heating is the same in every arm that still has
its luminosity.

**Four of the six arms leave the pile-up where it is, at the same amplitude, at the same
time.**  `c1`, `c2`, `d` and `f` change the top boundary, the top sponge, the taper window
and the radial resolution respectively, and the 0.97 R density at 1.5 turnovers is 2.38,
1.85, 2.38, 2.39 against the baseline's 2.36.  The two arms that do change it -- `a_mlt0`
and `e_handover` -- do so by taking the luminosity or its carrier away, so the shell is
STARVED instead of loaded; that is a different failure, not a cure.

**`f_nx144` is the one that settles the character of the defect.**  A truncation error
shrinks with the grid.  At `nx1 = 144` -- 1.5x the cells, the same stretch map, the
smallest cell 3.98e8 instead of 5.97e8 cm -- the pile-up is 2.39 against 2.36, `L_out/L`
is 0.960 against 0.970, and `dt` is SMALLER throughout.  It does not converge away, which
is what a FORMULATION error looks like: the grey spherical two-stream transports `J = A I`
with source `A B` (`he4_presn_cs.athinput`, the `rad_tau_lo` block), whose first moment
carries `1 - H_T/(2r)` of the exact flux at every resolution.

### 5.1 ARM (b): WITH THE WHOLE-COLUMN BLEND THERE IS NO CONDUCTION WALL TO INJECT INTO

`rt_bottom_flux = false` with `rad_flux_inner = 1.305278e15` was supposed to inject L
through the conduction wall instead of the two-stream's bottom face, as
`tests_3d/handover` did.  It injects **nothing**: under `rad_tau_lo/hi = 1e5/1e6` the blend
weight is `w_blend = 0` at every face (section 2.1), the radial conduction operator is
INERT, and `hydro/rad_flux_inner` is that operator's boundary flux.  The star is left with
no luminosity source at all: `L_out/L = 0.094` at half a turnover, the 0.97 R shell at 0.04
of its initial density and `T` at 0.38 of its initial value by one turnover.  **The arm is
not a bisection of the pile-up, it is a demonstration that `rt_bottom_flux = true` is
MANDATORY under the whole-column blend** -- and, conversely, that the `tests_3d/handover`
arms B/BH/D/DH were injecting L through a live operator only because they ran 20/300.

### 5.2 TWO ANOMALIES FOUND ON THE WAY, BOTH OPEN

1. **`c1_topwall` loses 34 % of the domain mass with BOTH radial boundaries walled and an
   EMPTY event log.**  `M/M0` is 1.0002 at 1.5 turnovers, then 0.915 (2.0), 0.7225 (2.5),
   0.6644 (3.0), while `he4.log` has not a single non-zero counter -- no `dfloor`, no
   `efloor`, no `vceil`, no `fofc`, no `eos_fail` -- and the baseline with an OPEN top over
   the same interval loses only 1.6 %.  A walled star cannot lose mass through a boundary
   and no floor fired, so this is a conservation violation of the interior scheme once the
   envelope goes sonic against the wall.  The first suspects are the dynamic well-balanced
   flux (`wellbalance_dynamic`, `wb_x1`, `wb_cache_every = 0`) at a wall face, and the
   deferred cubed-sphere floor path (`defer_cons_floors`, whose corrections are applied in
   `GnomonicEquiangleRaiseVel` and may not all be counted).  It does NOT affect the
   pile-up verdict -- `c1` reproduces the pile-up exactly over 0 to 1.5 turnovers, where
   its mass is conserved to 2e-4 -- but it needs its own round.
2. **The top sponge is holding `dt` up, not down.**  `c2_nosponge` has `dt = 0.561 s` at
   half a turnover where the baseline has 12.4 s, and yet the BEST emergent luminosity of
   any arm (1.003 and 0.993 L).  So `vdamp_top_tau = 0.1` is currently buying a factor 20
   in timestep by damping the very top cells -- i.e. by hiding the 3-7 s radiative sink of
   section 2.2 rather than by curing it.

## 6. THE VERDICT

1. **NUMERICAL, and the operator is named: the grey spherical two-stream's radial flux.**
   Not the taper (arm d moves the window and the flux profile does not change in the fifth
   significant figure), not the MLT closure (it carries zero above 0.94 R at t = 0, and
   switching it off replaces the pile with a drain), not the top boundary (a wall changes
   nothing over the first 1.5 turnovers), not the top sponge, and not the radial grid
   (1.5x resolution reproduces it to 1 %).  The two measurable defects are
   * `L_out/L = 0.956-0.981` once relaxed -- the documented `1 - H_T/(2 r)` spherical-form
     error, applied to all 96 cells because `rad_tau_lo/hi = 1e5/1e6` gives the two-stream
     the whole column -- against a domain thermal time of **0.462 turnover**, so the
     envelope banks 4-9 % of its own energy per turnover and `x_gate5` measures +20 % by
     2.5 turnovers;
   * a flux divergence in the last five cells (`F_2s/F_req` 1.007 -> 1.28) whose radiative
     cooling time is **3.4-6.5 s against `dt` = 12.9 s**, which empties the 0.99-1.01 R
     cells in the first half turnover and drops their mass onto 0.98 R.
   The PHYSICAL part is real but secondary: this envelope IS super-Eddington
   (`kappa F/(c g) = 1.23` at the MLT peak) and a 1-D column cannot convect, so once the
   banked energy pushes `Gamma(0.94 R)` through 1.0 -- measured, at exactly one turnover --
   it MUST lift.  The lifting is the physics; what makes it happen at one turnover instead
   of never is the 2-4 % of L the operator loses.
2. **No input change fixes this**, and that is the finding.  Of the five switches that
   were plausible, four do nothing and the fifth (`rad_tau 20/300`) makes it worse for the
   third time in a row.
3. **RECOMMENDED, in order.**
   * **Fix the spherical form of the two-stream**, which is the only change that can make
     `L_out/L = 1` and is the item `he4_presn_cs.athinput` has been carrying as a comment
     since the handover round: transport `I` with the `2 mu I/r` term treated as the
     geometry it is (or carry the first moment with an Eddington factor) so that the
     scheme's flux is the exact one rather than `1 - H_T/(2r)` of it.  A 1-D grey static
     atmosphere with an analytic answer is the test, and `L_out/L` off `rt_surface.bin` at
     0.1 turnover is the number to watch: it must be 1.000, not 0.98.
   * **Build the ic in the SOLVER's closure, not in an Eddington one.**  The 1.278 -> 0.981
     step in `L_out/L` within the first dump interval is the ic and the solver disagreeing
     by 28 % at the top face, and it is what kicks the 3.4 s cell.  One zero-cycle probe
     plus a secant pass on `T(r)` against the run's own `rt_surface`/`rt_profile` would
     close it, as `make_ic_he.py` already does for the box's pressure.
   * **Until either is done, do not raise `nx1` and do not change the taper**: both were
     measured here and neither moves the answer.  `rt_bottom_flux = true` must stay true
     as long as the blend is whole-column (5.1), and `vdamp_top_tau` must stay on (5.2),
     for the same bad reason.
   * The two open anomalies of 5.2 -- the walled-star mass loss and the sponge's hold on
     `dt` -- are each worth one short round, and the first is a correctness question about
     the interior scheme rather than about this star.

## 7. FILES

* `README.md` -- this file.
* `r9_arm.sh` -- one bisection arm; `base9/ a_mlt0/ b_botflux/ c1_topwall/ c2_nosponge/
  d_taperup/ e_handover/ f_nx144/` keep the `.hst`, the event log, `column_*.txt`,
  `mltfaces_*.txt`, `rt_profile.bin` and `rt_surface.bin`.  The `.log` beside each is the
  job's stdout.
* `r9_slotchk.sh` -- the 12 min `apudev` job on the fully fixed binary; `slotchk/`.
* `prof9.py` -- the shell table of section 2.3 (rho, T, v1, Mdot, tau, `F_2s/F_req`, the
  taper `w`, `Gamma`, `eint` at seven radii per 0.25 turnover) and the 0.97 R shell
  budget.  `--raw` corrects BOTH pre-fix slots, `--traw` only slot 5 (files from
  `athena_v5`).
* `budget9.py` -- the t = 0 face budget and shell divergence of sections 2.1 and 2.2, off
  the run's own `mltfaces` dump, plus the measured `d rho/dt` and `d eint/dt`.
* `lout9.py` -- `L_out/L` off `rt_surface.bin` (section 3).
* `arms9.py` -- the coarse arm table.
* `make_ic_sph_d.py` -- `bench/hestar_presn/make_ic_sph.py` with the taper window read off
  tau 1 / 0.1; writes `ic_he4_presn_sph_tau1.txt` HERE and does not touch the committed
  ic.  Its window: `eos_rad_rho_hi = 9.0918e-10`, `rho_lo = 2.7774e-10`,
  `t_hi = 5.18187e4`, `t_lo = 4.26635e4`.
* `build.sh`, `build_gpu.log`, `build_gpu2.log` -- the two builds.  `athena_v5` has only
  the slot-6 fix (a stale object survived the first build's kill; the arms were launched
  with it, which is why `prof9.py` has `--traw`); **`athena_v6` is the binary to keep**.
  Neither is committed; rebuild from this commit in `build_gpu_rg`.

## 8. Coordinator's caveat (added 2026-09-18 01:00)

Section 6 attributes the 2-4 % emergent-luminosity shortfall to the `1 - H_T/(2r)` spherical-form
error. That error was removed in commit 8bca3dfa (tests_r4): the thick unit test is exact to 1e-4
over the interior. The measured shortfall sits in the last five cells (dtau_cell < 0.3) and at the
top face (`F_2s/F_req` = 1.28 at t = 0), i.e. in the optically thin layers and the top boundary
condition of the spherical sweep (rt_top_re off: the ghost mirrors the top cell; the thick test
already showed 1e-3..1e-2 residuals at face 122). The recommended test stands: a 1-D grey
atmosphere with an analytic answer, gating on L_out/L = 1.000 off rt_surface.bin at 0.1 turnover,
but the suspect is the top-face / thin-layer treatment, and the ic-vs-solver mismatch at the top
(1.278 -> 0.981 within the first dump) is the second item.
