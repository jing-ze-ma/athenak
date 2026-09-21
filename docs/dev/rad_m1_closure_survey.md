# M1 Eddington-factor closures for grey photon RHD of a stellar envelope

Literature survey, 2026-09-21. Question: which `chi(f)` should `src/rad_m1/rad_m1_closure.hpp`
offer next to Levermore (1984), for a stellar envelope/atmosphere (thick CZ below, a
plane-parallel-like diffuse photosphere, thin above, `f ~ 0.3-0.7` at the surface, no beams),
solved implicitly with a Picard-lagged closure?

Everything marked **[computed]** I evaluated myself from the published formulas (script kept in the
session scratchpad, not in the repo). Everything marked **[unverified]** I could not confirm from a
source I actually read.

---

## 1. The closures, with formulas

`f = |F|/(cE)` (= `H/J`), `chi = P/E` (= `K/J`), `P_ij = E[(1-chi)/2 d_ij + (3chi-1)/2 n_i n_j]`.

| Closure | `chi(f)` | Reference |
| --- | --- | --- |
| Eddington | `1/3` | classical |
| **Levermore 1984** (boosted isotropic) | `(3 + 4f^2) / (5 + 2 sqrt(4 - 3f^2))` | Levermore, *JQSRT* **31**, 149 (1984), "Relating Eddington factors to flux limiters" |
| **Minerbo 1978** (classical max. entropy) | `1/3 + (2 f^2/15)(3 - f + 3 f^2)` | Minerbo, *JQSRT* **20**, 541 (1978) |
| **Kershaw 1976** | `(1 + 2 f^2)/3` | Kershaw, LLNL report UCRL-78378 (1976) |
| **Wilson** (Wilson et al. 1975 / Bruenn) | `1/3 - f/3 + f^2` | quoted as "Wilson" in Murchikova+2017 eq. for Wilson et al. (1975) |
| **Janka 1991/1992** (Monte-Carlo fit) | `(1/3)[1 + a f^m + (2-a) f^n]`, sets `{a=0.5, m=1.3064, n=4.1342}` (J1) and `{a=1, m=1.345, n=5.1717}` (J2) | Janka, PhD thesis / A&A 1991-92; parameters as tabulated by Murchikova+2017 |
| **Cernohorsky & Bludman 1994** (Fermi-Dirac ME, "MEFD") | `1/3 + (2/3)(1-e)(1-2e) X(f/(1-e))`, with the polynomial `X(x) = x^2 (3 - x + 3x^2)/5` and `e` the mean occupation | Cernohorsky & Bludman, *ApJ* **433**, 250 (1994) |
| **Levermore-Pomraning 1981** (flux-limiter-derived) | implicit: `lambda(R) = (2+R)/(6+3R+R^2)` (rational fit) or `lambda = (1/R)(coth R - 1/R)` (exact); then `f = lambda R` and `chi = lambda + lambda^2 R^2 = lambda + f^2` | Levermore & Pomraning, *ApJ* **248**, 321 (1981); the `chi = lambda + lambda^2 R^2` identification is Levermore (1984) |

Notes.

* MEFD reduces **exactly** to Minerbo as `e -> 0`. For **photons** the relevant max-entropy family is
  Bose-Einstein, not Fermi-Dirac; the Bose-Einstein ME closure is the Minerbo form in the
  classical (`e << 1`) limit, and I could not find a published closed-form BE Eddington factor for
  finite occupation that is used in production photon codes — **[unverified]**. For a stellar
  envelope the photon occupation is not small at depth, but there `f -> 0` and every closure
  collapses to 1/3, so this is moot for us.
* LP81 is a *flux limiter*, not a closure; its induced `chi(f)` is what Levermore 1984 derived, and
  it is the worst of the set on the atmosphere test (sec. 3).

---

## 2. Question 1 — the grey plane-parallel (Milne/Hopf) atmosphere

Exact, from the moment equations of the conservative grey half-space (`dH/dtau = 0`,
`dK/dtau = H`, `J = 3H(tau + q(tau))` with Hopf's `q`):

* `K(tau) = H (tau + q(inf))` exactly, so at the surface `K(0) = H q(inf)`.
* Hopf's exact endpoint values: `q(0) = 1/sqrt(3) = 0.5773503` and `q(inf) = 0.7104460895`.

Therefore **[computed]**

```
f(0)   = H/J  = 1/(3 q(0)) = 1/sqrt(3) = 0.577350
chi(0) = K/J  = q(inf)/(3 q(0)) = q(inf)/sqrt(3) = 0.410176
```

(This is the familiar "0.4101 at the surface, 1/3 at depth" statement. I derived it from `q(0)` and
`q(inf)`; I did not find a page that prints 0.410176 itself, so treat the *derivation* as mine and
the two Hopf constants as the cited inputs.)

What each closure predicts at `f = 0.577350` **[computed]**:

| closure | `chi(f0)` | error vs 0.4102 | `f` at which it gives 0.4102 |
| --- | --- | --- | --- |
| Eddington `1/3` | 0.3333 | **-0.0768** | — |
| Janka 1 (`a=0.5`) | 0.4663 | +0.0561 | 0.449 |
| Wilson | 0.4742 | +0.0640 | 0.490 |
| **Minerbo 1978** | **0.4855** | **+0.0753** | 0.430 |
| **Levermore 1984** | **0.5120** | **+0.1018** | 0.386 |
| Janka 2 (`a=1`) | 0.5120 | +0.1018 | 0.332 |
| Kershaw 1976 | 0.5556 | +0.1454 | 0.340 |
| LP81 (rational `lambda`) | 0.5586 | +0.1484 | 0.346 |
| LP81 (exact `lambda`) | 0.5904 | +0.1803 | — |

**Every `f`-only closure overshoots**, because the true surface field is a quasi-isotropic outgoing
half-space distribution (`I > 0` only for `mu > 0`, but broad in `mu`), whereas an M1-type closure
with `f = 0.577` reconstructs something much more forward-peaked. Ranking by `|error|` at the
atmosphere point: Janka 1 < Wilson < **Minerbo** < Levermore = Janka 2 < Kershaw < LP81. Plain
Eddington `1/3` errs by -0.077, i.e. *less* than Levermore's +0.102 and comparable to Minerbo's
+0.075 — worth knowing: for this one problem, M1-Levermore is worse than pure Eddington.

`chi` over the band that matters to us **[computed]**:

| `f` | 0.3 | 0.5 | 0.7 | 0.9 |
| --- | --- | --- | --- | --- |
| Levermore | 0.3791 | 0.4648 | 0.6063 | 0.8313 |
| Minerbo | 0.3690 | 0.4417 | 0.5796 | 0.8226 |
| Kershaw | 0.3933 | 0.5000 | 0.6600 | 0.8733 |
| Wilson | 0.3233 | 0.4167 | 0.5900 | 0.8433 |
| Janka 1 | 0.3714 | 0.4292 | 0.5524 | 0.8020 |

Minerbo sits ~2-5 % below Levermore through the whole `f = 0.3-0.7` window; Janka 1 a further
1-3 % below that. None of them is close to 0.41 at `f = 0.58`: **the residual error is structural to
any closure that is a function of `f` alone**, and only a half-space/half-moment or a genuine VET
method removes it (sec. 4).

---

## 3. Question 2 — pathologies, realizability, smoothness

**Two-beam / crossing beams.** The canonical M1 failure: two equal counter-propagating beams give
`F = 0`, hence `f = 0`, hence `chi = 1/3` for *every* `f`-only closure, whereas the true `chi = 1`.
M1 then makes the beams collide and produce a spurious stagnation/shock instead of passing through.
Documented for photons (Rosdahl & Teyssier 2015, RAMSES-RT: "fails in-between multiple sources,
creating spurious sources of perpendicular radiation"), for neutrinos (Foucart et al. 2018,
*MNRAS* **475**, 4186 — Monte-Carlo closure comparison; M1 "very inaccurate in the low-density polar
regions", spurious neutrino shocks on the axis), and for disks (Melon Fuksman et al. 2025, below:
M1 midplane temperature 44 % too high at one group, 21 % at 22 groups, because vertical inward and
outward fluxes merge). **This is exactly our geometry**: in a diffuse photosphere and in the thin
region above a corrugated convective surface, up- and down-going radiation coexist.

**Shadows.** M1 casts a shadow; Eddington/FLD does not. Menon et al. 2022 (VETTAM, *MNRAS* **512**,
401) state that neither Eddington nor M1 casts a *qualitatively correct* shadow, whereas VET does.
Anninos & Fragile 2020 (*ApJ* **900**, 71) cite shadow preservation as M1's main advantage over
Eddington. Shadow quality is not a discriminator among `chi(f)` choices — all M1 forms shadow — so
it is a weak test for our purpose.

**Realizability.** For an intensity `I >= 0` on `mu in [-1,1]`, `f^2 <= chi <= 1`. **[computed]** all
of Levermore, Minerbo, Kershaw and Wilson satisfy `chi >= f^2` and `chi <= 1` on `0 <= f <= 1`, so
on the *scalar* bounds none of them is violated. The stronger statement usually made is that
Levermore/M1 and Minerbo/ME are realizable *by construction* (they come from an actual non-negative
distribution) whereas Kershaw and Wilson are interpolations that merely happen to satisfy the
bounds. I could **not** find a source stating that Wilson is non-realizable — **[unverified]**;
what *is* true and checkable is that Wilson has `dchi/df = -1/3` at `f = 0` and `chi < 1/3` for
`f < 1/3` (e.g. `chi(0.3) = 0.3233`), i.e. it dips *below* the isotropic value, which is
unphysical-looking and will fight a Picard iteration near the diffusion limit.

**Hyperbolicity / characteristic speeds.** Levermore 1984 gives the M1 characteristic speeds; the
system is hyperbolic with speeds `<= c` for the ME/Levermore closure. Sec. 2 of
arXiv:2603.22400 ("Second order closures for the radiative transfer equation: some are unstable",
2026) confirms M1 (Levermore) and the maximum-entropy closure are stable hyperbolic two-moment
formulations at first order (their stability condition is `(A f)' <= 1` with `A = chi`), and shows
that *second-order* local closures and OTVET extensions are unstable. Kershaw's linear-in-`f^2`
form has the smallest `dchi/df` at `f -> 1` (**[computed]** 1.33 vs 1.99 for Levermore/Minerbo),
so its free-streaming characteristic speed is least sharp — one reason it does badly in transparent
regions. Fermionic Kershaw-type closures with the realizability bounds built in exist
(Schneider & Alldredge, *ZAMP* 68 (2017), "Kershaw-type transport equations for fermionic
radiation") — a tunable convex family — but I did not find a photon/atmosphere application
**[unverified]**.

**Smoothness — this is what matters for our lagged closure.** **[computed]**

| closure | `f -> 0` behaviour | `dchi/df` at `f -> 1` |
| --- | --- | --- |
| Levermore | `1/3 + f^2/2 + O(f^4)`, analytic | 1.99 |
| Minerbo | `1/3 + (2/5) f^2 - (2/15) f^3 + (2/5) f^4`, **polynomial, exactly analytic** | 1.99 |
| Kershaw | `1/3 + (2/3) f^2`, polynomial | 1.33 |
| Wilson | `1/3 - f/3 + f^2`, slope `-1/3` at 0 | 1.66 |
| Janka 1/2 | `1/3 + (a/3) f^m`, `m ~ 1.3`: **fractional power, `chi''` diverges at `f = 0`** | 2.28 / 2.16 |
| LP81 (exact) | non-analytic (`coth`), requires a root find per cell | — |

For a Picard-lagged closure the ordering is clear: **Minerbo and Kershaw are polynomials** (no
`sqrt`, no branch, `chi` and `dchi/df` bounded and smooth everywhere including `f = 0` and `f = 1`),
Levermore has a `sqrt(4-3f^2)` that is smooth but needs the `f <= 1` guard we already have, Wilson
has the wrong sign of `dchi/df` at small `f`, and Janka's fractional powers are both expensive
(`pow`) and `C^1`-but-not-`C^2` at `f = 0`, which is where most of our cells live. Janka is
therefore ruled out for us despite its good value at the atmosphere point.

---

## 4. Question 3 — closures designed for or tested in atmospheres / half-spaces

* **Half-moment (HM) closure — Melon Fuksman, Flock, Klahr, Mattia & Muley 2025**
  (arXiv:2504.13999, "Multidimensional half-moment multigroup radiative transfer"). Carries **two**
  M1-like systems, one per hemisphere along a preferred direction (the disk vertical / our radial
  direction), so up- and down-going radiation never merge. Against Monte Carlo in a protoplanetary
  disk: M1 midplane temperature +44 % (1 group) / +21 % (22 groups); HM +6 % (22 groups), and with
  only 3 groups HM's worst departure is 8 % vs M1's 23 %. This is the single most directly relevant
  result to a plane-parallel photosphere, and it is a *geometry* fix, not a `chi(f)` fix. Cost is
  roughly 2x the moments.
* **Levermore-Pomraning / FLD family** was built for exactly this problem and is the *worst* of the
  set at the Hopf point (sec. 2). Not a route.
* **Optical-depth-dependent interpolation between Eddington and M1** — I looked for a published,
  tested "blend `chi` towards 1/3 as a function of local `tau`" closure and did **not** find one as
  a named scheme **[unverified]**. Our module already does something of this kind in the *flux
  limiter* (the `F_diff` / `alpha` blend in `rad_m1_closure.hpp` deliberately uses `1/3`), which is
  the right place for it; an ad-hoc `chi` blend would break realizability and the hyperbolic
  structure.
* **Cernohorsky & Bludman 1994 MEFD** uses extra information (the occupation `e`). Murchikova+2017
  found ME and MEFD "yield almost identical results" under PNS conditions and prefer ME for
  simplicity. For photons the analogue would be Bose-Einstein ME; no production photon code I found
  uses it **[unverified]**.
* **Machine-learned closures** exist and beat analytic ones in the trained regime — Harada et al.
  (Deep learning of the Eddington tensor in CCSN, arXiv:2104.13039), a LightGBM extended closure
  (*ApJ* 2025, arXiv:2409.02719), and "A hyperbolic neural closure for M1 radiation transfer"
  (arXiv:2607.10364), which exists precisely because unconstrained learned closures produce
  non-real characteristic speeds. Not appropriate for us (table/inference cost, no training data,
  and the task brief asks for formula-only options).
* **VET (Jiang, Stone & Davis 2012, *ApJS* **199**, 14; Menon et al. 2022 VETTAM; Jiang's later
  Athena++ VET)** is the reference answer if one wants the atmosphere right: the Eddington tensor
  comes from a short-characteristics formal solution, so the 0.4102 surface value is reproduced by
  construction. Non-local, needs a global angular solve — out of scope for this module.
* **Photon-RHD codes that chose Levermore anyway**: Skinner & Ostriker 2013 (*ApJS* **206**, 21,
  Athena two-moment module), QUOKKA (Wibking & Krumholz 2022, *MNRAS* **512**, 1430, Levermore by
  default), RAMSES-RT (Rosdahl & Teyssier 2015), Melon Fuksman & Mignone 2019 (*ApJS* **242**, 20,
  PLUTO), Anninos & Fragile 2020. None of them reports a grey-atmosphere Hopf test. So there is no
  published evidence that any *other* `chi(f)` is preferred for photon envelopes — the published
  comparisons that do rank closures (Murchikova+2017; Richers 2020, *PRD* **102**, 083017) are
  neutrino problems.
* **Murchikova, Abdikamalov & Urbatsch 2017** (*MNRAS* **469**, 1725, arXiv:1701.07027) vs Monte
  Carlo around PNS models: **Minerbo "on average yields relatively accurate results in the broadest
  set of cases"** and is never the worst; no closure wins everywhere; **Wilson and Levermore are
  best in transparent regions, Kershaw in opaque zones**. **Richers 2020**: of 8 closures, MEFD,
  Levermore and Janka 2 are "all reasonable", but the assumptions underlying *all* analytic closures
  are strongly violated at thick/thin interfaces — which is our photosphere.
* **Olson, Auer & Hall 2000** (*JQSRT* **64**, 619) compare diffusion, P1, several FLD theories and
  VEF closures on two time-dependent test problems (the Marshak-type problems that became the
  standard "Olson-Auer-Hall" benchmarks). I could not fetch the full text (OSTI/publisher blocked),
  so I have **not** verified which closure they rank best — **[unverified]**.
* **Smit, van den Horn & Bludman 2000** (*A&A* **356**, 559, "Closure in flux-limited neutrino
  diffusion and two-moment transport"). Cited widely as the reference that the maximum-entropy
  closure family is the right one to use in the two-moment system and on the relation between flux
  limiters and Eddington factors. I could not fetch the text; its detailed conclusions are
  **[unverified]** here.

---

## 5. Question 4 — recommendation

**Implement two, both one-liners, selected by an input flag alongside `levermore` and `eddington`:**

1. **Minerbo 1978 / maximum entropy** — `chi = 1/3 + (2 f^2/15)(3 - f + 3 f^2)`.
   Rationale: it is the closure with the best published broad-spectrum record against Monte Carlo
   (Murchikova+2017), it is a pure polynomial (cheapest of all, no `sqrt`, `C^inf` on `[0,1]`,
   analytic `dchi/df = (2/15)(6f - 3f^2 + 12 f^3)` if we ever want a Newton closure instead of
   Picard), it is realizable by construction, and **[computed]** it is 25 % closer than Levermore to
   the exact Hopf surface value (+0.075 vs +0.102) and sits 2-5 % below Levermore through the whole
   `f = 0.3-0.7` band we care about. Lower `chi` at fixed `f` means less spurious radiation
   *pressure anisotropy* in the photosphere, which is the quantity that couples to the gas in our
   envelope runs.
2. **Kershaw 1976** — `chi = (1 + 2 f^2)/3`.
   Not because it is accurate (it is the worst analytic one at the Hopf point) but because it is the
   *cheapest possible* smooth closure with the correct limits, it is the best in opaque zones per
   Murchikova+2017, it has the smallest `dchi/df` at `f -> 1` so it is the most benign for the
   implicit solve, and it brackets Levermore/Minerbo from above. Having Eddington (1/3), Minerbo,
   Levermore and Kershaw gives a four-point spread that turns "closure sensitivity" into a
   measurable number for any result we publish.

Do **not** implement: Wilson (negative `dchi/df` at small `f`), Janka (fractional powers, `pow` cost,
`C^2`-singular at `f = 0`), LP81 (needs a root find, worst at the Hopf point), MEFD (fermionic;
degenerates to Minerbo for photons), ML closures.

**Discriminating test — the steady grey atmosphere gate is the one that matters.** Run it with all
four closures and report, at the top of the atmosphere, the pair `(f, chi)` and the emergent
`T(tau)` / Hopf function against the exact `q(tau)`:

* the gate's closed form already gives `q(0) = 1/sqrt(3)`, `q(inf) = 0.71045`, so the expected
  surface point is `(0.5774, 0.4102)`;
* prediction **[computed]**: Eddington will land at `chi = 1/3` (`-0.077`), Minerbo at `0.486`
  (`+0.075`), Levermore at `0.512` (`+0.102`), Kershaw at `0.556` (`+0.145`), and the surface
  temperature / limb behaviour should order the same way. If the measured spread is smaller than
  this, the gate is not actually sampling `f ~ 0.58` and the test needs a thinner top.
* The **shadow test will not discriminate** (all four shadow; Eddington does not) and the
  **Marshak wave barely will** (it lives at `f -> 0` where all four agree to `O(f^2)`, differing only
  in the `f^2` coefficient: 1/2 Levermore, 2/5 Minerbo, 2/3 Kershaw). **Radiative shocks** probe the
  intermediate `f` and are a useful secondary check, mainly on whether the Picard iteration count
  changes.
* Free bonus check: a 1-D two-beam (crossing-beam) setup will fail identically for all four
  (`f = 0 -> chi = 1/3`). If that failure mode turns out to matter for the seeded 2-D He slab
  (up-going and down-going flux coexisting above the photosphere), the fix is **not** a different
  `chi(f)` but the half-moment formulation of Melon Fuksman et al. 2025 — worth keeping on the list
  as a later, larger change.
