# tests_ck_sweep_form — the spherical face coupling solved EXACTLY, in two passes

`problem/ck_sweep_form` in `src/utils/two_stream_rt.hpp`, read in
`src/pgen/deep_hot_jupiter_rt.cpp`.  Follow-on to `tests_ck_sph` (the spherical form
itself) and `tests_ck_sweep_cost` (what it costs).

## 1. What the four passes were for, and why two are enough

Under `ck_spherical` the up and down rays are coupled at every face by

    c = beta (d_above - u_below),   beta = (A_above - A_below)/(A_above + A_below),
    u_above = u_below + c,          d_below = d_above + c,

so a one-directional down-sweep needs `u_below`, which it does not have.  The shipped
kernel estimates it with TWO PROBE PASSES (a plane-parallel down probe `d0`, then an up
probe carrying the mixing against `d0`) and only then runs the real down and up sweeps:
**four column passes per chain**, and a `c` that is right only to `O(beta^2)`.

The structure the probes were working around is simpler than it looks.  Every half layer
in the column **transmits with no reflection** — pure absorption, no scattering — and the
ONLY reflection in the whole column is the face itself, of amplitude `beta`:

    u_above = (1-beta) u_below + beta d_above,
    d_below = (1+beta) d_above - beta u_below.

That is a stack of transmitting layers separated by reflecting interfaces, i.e. the
textbook adding / invariant-imbedding problem, whose exact solution is **two passes**:
carry the linear relation between the two rays up from the cut face, impose the top
boundary datum, come back down.  `ck_sweep_form` does exactly that, twice over.

### 1.1 `ck_sweep_form = 1` — `tm`, the transfer-matrix / adding form, in (u, d)

Carry, at each face, `u_below = R d_below + Sc` with `R` the reflectance of everything
below the face and `Sc` what it emits.  A half layer of transmission `T = 1 - e0`, whose
emission along an up crossing is `P` and along a down crossing is `Q`, maps

    R  <-  T^2 R,        Sc  <-  T (R Q + Sc) + P,

and the face is the reflection-addition Moebius map

    R  <-  (R + beta)/(1 + R beta),     Sc  <-  (1 - beta) Sc/(1 + R beta).

**Conditioning**: no raw product of exponentials appears anywhere — that is the whole
reason for carrying reflectances rather than 2x2 transfer matrices.  A thick layer sends
`T^2 R` to zero harmlessly, where a transfer matrix would need `1/T = e^{+tau}` and
overflow on the first `tau ~ 700` cell (this column reaches `tau ~ 1e4` at the cut).
`|R| < 1` and `|beta| < 1`, so `1 + R beta` is bounded away from zero.

Cost: **no divide in the layer**, one divide at the face.

### 1.2 `ck_sweep_form = 2` — `sd`, in S = (u+d)/2 and D = (u-d)/2

The same statement written in the variables in which the face conditions are
"`S` passes unchanged, `A D` passes unchanged", i.e. **diagonal, no cross-coupling**.
Writing the carried relation as `D = -G S + H`, which is `R`, `Sc` re-expressed through

    G = (1-R)/(1+R),     H = Sc/(1+R),

the face becomes a pure RESCALING,

    G <- rho G,   H <- rho H,    rho = A_below/A_above = (1-beta)/(1+beta),

two multiplies and no divide, which IS the statement that (S, D) diagonalise the face.
The cost moves into the layer, which now needs one divide per half layer:

    den = (1+G) + T^2 (1-G)
    G  <-  [(1-T^2) + G(1+T^2)]/den
    H  <-  [T((1-G) Q + 2H) + P(1+G)]/den

Two half layers per cell against one face per cell, so **sd trades one divide for two**:
it is the more elegant form and the more expensive one.  Section 4 measures that.

### 1.3 Why one cannot simply "sweep S and D"

Inside a layer the two-stream system is `dS/dt = -D`, `dD/dt = -(S - B)`, i.e.
`S'' = S - B`, whose solutions are `e^{+t}` as well as `e^{-t}`.  A forward recurrence in
(S, D) amplifies like `e^{+tau}` and dies on the first optically thick cell; `u` and `d`
ARE the characteristics, and what (S, D) buys is the FACE, not the layer.  The Riccati
variable `G` above is the stable way to carry the same statement — and it is the one that
makes the face free.

### 1.4 The boundary data

* **Bottom (the correlated-k cut face).**  `u_below` is given: the cut cell's Planck
  function plus the internal flux, in the cut FACE's frame.  So the relation starts at
  `R = 0`, `Sc = B_cut + Iint` — in (S, D), `G = 1`, `H = B_cut + Iint`, i.e. `S + D = u`.
  `G = 1` is a perfectly absorbing bottom; a face is the only thing that can make `G < 1`.
* **Top.**  `d_above` is given: the unresolved hydrostatic ghost column's emission
  `(1 - e^{-dtau}) B_ghost`, in the top FACE's own frame.  Pass 2 opens by combining it
  with the accumulated relation — one divide — and everything else follows downward.
* **Both are unchanged data**; nothing about the boundary conditions had to be
  re-derived, only re-expressed.

### 1.5 The deposit still telescopes — with nothing left over

Both forms enforce `A D` continuity the same way the four-pass form does, by
CONSTRUCTION rather than by inference: at each face the flux is reported as
`A_below (u_b - d_b)/A_face` (per unit area at its own face) and the up ray is handed
`u_above = d_above + (u_b - d_b) A_below/A_above`.  Cell `i` then receives, from the two
rays' four half-layer crossings,

    Src_i dx_i / w = [d_b(i+1) - d_a(i)] + [u_a(i) - u_b(i+1)]
                   = 2 D_a(i) - 2 D_b(i+1)
                   = -(A_{i+1} F_{i+1} - A_i F_i)/ACC_i,

which is the required `-(A_t F_t - A_b F_b)/V` exactly.  In pass 2 the down ray is
STEPPED and the up ray is RECOVERED from the stored relation, so the up ray's exit and
the `u_below` the face above reported agree only to round-off; the sweep adds that
difference back into the same cell (one subtract, one add per cell and chain), which is
what makes the per-cell identity hold with nothing left over.

### 1.6 No new private memory

The probe form needs `I_down` (the whole down column) and `Cmx` (the face mixing).  The
two-pass form needs NEITHER — pass 2 has the down intensity and the up intensity at the
same face in the same iteration — so those two `NC x NN` columns are borrowed to hold the
Riccati pair instead.  The thread's private segment is the size it was.

## 2. What was touched

| file | change |
| --- | --- |
| `src/utils/two_stream_rt.hpp` | the flag and its note; a startup guard; `const int ckform_`; a `cofs` lambda next to `step` (the same coefficient triple, handed back instead of applied); the two-pass block inside `rt_chain_ck` under `if constexpr (FRM != 0)`, which returns before the four-pass code; a `FRM` tag on the kernel dispatch |
| `src/pgen/deep_hot_jupiter_rt.cpp` | `problem/ck_sweep_form` (0/1/2) and its range check |

The four-pass code is not touched by one character; the new form is a separate branch
reached only when the flag is set, and the flag is a COMPILE-TIME tag, so the default
instantiation is literally the code it was (gate a).

**Collision, declared.**  `launch_ck_chain` / `launch_ck_cache` are part of an
uncommitted `ck_sweep_cache` change by another agent; they gained one template tag each
(`auto frm_tag`, and the three call sites that pass it through) plus a new
`launch_ck_form` in between.  Mechanical, no behaviour of theirs changed.

**Refused at startup** (rather than silently changed): `ck_spherical = false` (with no
face coupling the plane-parallel sweep is already the exact two-pass form),
`ck_sweep_cache != 2` (the new forms are instantiated at the column cache only, to keep
this kernel's instantiation count at 8 per radial tier instead of 12), and `ck_implicit`
(its tridiagonal is the linearisation of the FOUR-PASS recurrence and has not been
re-derived for the Riccati one).

## 3. Gates

Serial CPU, `build_cksd` (`cmake -B build_cksd -D PROBLEM=deep_hot_jupiter_rt`,
Release), input `dhj_ck_sweep_form.athinput` = `inputs/tests/dhj_ck_spherical.athinput`
with `ck_spherical = true` and the new key added (a CLI override needs the key to exist).

### (a) The default form is BITWISE what HEAD runs — PASS

`a_head` = `446b7b30` (`git archive HEAD`, built clean in
`bench/cksd_0922/ref_head/build_cpu`), `a_new` = this tree at `ck_sweep_form = 0`.
20 cycles, `ck_spherical = true`, hst and a binary dump every cycle, plus a column dump:

| comparison | result |
| --- | --- |
| `dhj.hydro.hst` | byte-identical |
| column dump `col.txt` | byte-identical |
| 22 `bin/*.bin` data payloads | bitwise identical, 22 of 22 |

(The `.bin` FILES differ by 491 bytes — the embedded input text carries the extra
`ck_sweep_form` line.  The payload comparison is `bin_convert.read_binary` + `array_equal`
over every MeshBlock of every variable.)

This also re-confirms, from HEAD, that the uncommitted `ck_sweep_cache` change in the
same file is inert at its default.

### (b) tm and sd against the four-pass form — a REAL difference, and it is the probe

One-cycle column dumps (`problem/ck_dump_file`), three columns: day `mu0 = +0.922`
(`m=0,k=5`), twilight `+0.382` (`m=1,k=3`), night `-0.922` (`m=0,k=2`).  Max difference
normalised by the column maximum of the same quantity:

| pair | T | F_lw (net face flux) | Src_lw (deposit) | Q_sw (beam) |
| --- | --- | --- | --- | --- |
| tm vs sd | 0 | **3.5e-17** | **9.0e-13** | 0 |
| 4-pass vs tm | 0 | **2.33e-02** | **2.64e-02** | 0 |
| 4-pass vs sd | 0 | 2.33e-02 | 2.64e-02 | 0 |

identical on all three columns — as they must be, since the THERMAL two-stream does not
see `mu0` and the beam is untouched by either change (`Q_sw` differs by exactly zero).

After **20 cycles**, whole state, max difference normalised by the state maximum:

| pair | dens | eint | velx | vely | velz |
| --- | --- | --- | --- | --- | --- |
| tm vs sd | **0 (bitwise)** | 0 | 0 | 0 | 0 |
| 4-pass vs tm/sd | 1.4e-08 | 1.5e-07 | 6.9e-03 | 3.2e-04 | 5.9e-04 |

(the velocity is still near zero 20 cycles after a cold start, so its normalisation is
the least meaningful of the five.)

**So the answer to the question "is it round-off?" is NO, and the reason is that the two
are not the same linear problem.**  The four-pass form's `c` is built from a PROBE and is
accurate only to `O(beta^2)`; the two-pass forms solve the same coupling exactly.  Two
pieces of evidence that the new forms are the exact ones and the discrepancy is the
probe's lag:

1. **tm and sd agree to round-off.**  They are algebraically different routes — one
   carries a reflectance through a Moebius map at the face and two multiplies in the
   layer, the other carries an Eddington-like ratio through two multiplies at the face
   and a Moebius map in the layer — and they land on the same numbers to 3.5e-17 in the
   flux and 9e-13 in the deposit, and bitwise in the evolved state.  Two wrong answers do
   not do that.
2. **The discrepancy scales with the geometry it comes from.**  Shrinking `mesh/x1max`
   shrinks the correlated-k column's area ratio, and with it the whole coupling:

   | `x1max` | ck domain | `A_out/A_in` | 4-pass vs tm, F_lw | Src_lw |
   | --- | --- | --- | --- | --- |
   | 9.912e9 | empty (`icut > ie`) | — | **0** | 0 |
   | 1.0856e10 | 1 cell | — | **0** | 0 |
   | 1.3216e10 | r 1.108e10..1.322e10 | 1.422 | 6.5e-05 | 9.3e-05 |
   | 2.0556e10 (production) | r 1.107e10..2.056e10 | 3.446 | 2.33e-02 | 2.64e-02 |

   i.e. it vanishes identically where there is no face coupling to get wrong, and grows
   steeply with the area ratio — faster than `beta^2` per face because ~100 faces
   accumulate and `beta` itself grows outward.

### (c) The thermal gates of `tests_ck_sph`, with each form — PASS

`tests_ck_sph/budget2.py` on the day column (`mu0 = 0.922`, `A_out/A_in = 3.446`):

| form | per-cell flux identity `max|V dep - dPhi|/max|Phi|` | column budget `[sum V Src + (A_t F_t - A_b F_b)]/max|A F|` |
| --- | --- | --- |
| 0 (four-pass) | 2.808e-10 | 3.871e-09 |
| 1 (tm) | **2.822e-10** | **3.864e-09** |
| 2 (sd) | **2.822e-10** | **3.864e-09** |

i.e. all three conserve at the same level, which is the `RtF` round-off of the column
(`tests_ck_sph` quotes 2.05e-10 / 2.62e-09 for its own column).  The BEAM diagnostics of
the same script — absorbed power, `tau_slant = 1` radius, the power entering the top face
— are identical to all printed digits across the three forms, as they must be.

**The star-off `L(r)` flat gate.**  `g_relax_f{0,1,2}` = 4000-cycle relaxations with each
form; `rerun_night.sh` then restarts each from its last restart and dumps the NIGHT
column (`mu0 = -0.922`, no beam), and `tests_ck_sph/lprof.py` reads them.  Over the 52
faces the two-stream owns well above the radiative-conduction handover
(`r > 1.4e10`, `w_diff = 0`), the spread of `L = A F`:

| form | (max - min)/mean of L(r) |
| --- | --- |
| 0 (four-pass) | 3.073e-03 |
| 1 (tm) | 3.077e-03 |
| 2 (sd) | 3.077e-03 |

All three are flat at the same level, and none shows anything like the factor
`A_out/A_in = 2.47` by which the PLANE-PARALLEL form fails this gate
(`tests_ck_sph/README.md` gate 3).  Caveat, stated because it bounds the gate: 4000
cycles is a partial relaxation — `L(top)/L_int` is 992, not 1 — so this says the two-pass
forms conserve and dilute exactly as the four-pass one does, not that either is closer to
the true equilibrium.  Distinguishing THOSE would need a fully relaxed pair, which is a
much longer run than this directory holds.

## 4. GPU cost

`bench/cksd_0922`: HIP build of an rsync snapshot of this tree
(`bench/cksd_0922/build.sh`, recipe copied from `bench/cksw_0922/build.sh`), the restart
`bench/cs_mhd_prod3/rst/dhj.00567.rst` (rot 283, `ncycle = 4429831`), 300 cycles via
`nlim = 4430131`, 2 ranks / 2 GPUs on `apudev`, with
`bench/prof_0922/plain/{deep_hot_jupiter.athinput,submit.sh}` verbatim apart from `-J`,
the binary and the single `problem/ck_sweep_*` override (the two keys are appended to the
input's last `<problem>` block, because a CLI override needs the key to exist).

### 4.1 Registers and private memory, from the code-object metadata

`bench/cksd_0922/regs.sh`, i.e. `llvm-objcopy --dump-section=.hip_fatbin` on the pgen
object, `clang-offload-bundler --unbundle --targets=hipv4-amdgcn-amd-amdhsa--gfx942`,
`llvm-readelf --notes`.  Production tier: `NN = 136` (`n1 = 128`), `SPH = true`,
`BSP = false`.

| instantiation (NN, SPH, BSP, CCH, FRM) | vgpr | sgpr | agpr | vgpr spill | private segment |
| --- | --- | --- | --- | --- | --- |
| 136, 1, 0, 0, 0 — four-pass, cache off | 128 | 108 | 64 | 0 | 10048 B |
| 136, 1, 0, 1, 0 — four-pass, carry | 128 | 108 | 64 | 0 | 10192 B |
| 136, 1, 0, 2, 0 — four-pass + column cache | 128 | 108 | 64 | 0 | 27360 B |
| 136, 1, 0, 2, 1 — **tm** | 128 | 108 | 64 | 0 | **27136 B** |
| 136, 1, 0, 2, 2 — **sd** | 128 | 108 | 64 | 0 | **27232 B** |

**Registers do not discriminate.**  Every instantiation in the binary reports vgpr 128 /
sgpr 108 / agpr 64 with `max_flat_workgroup_size 1024` — the compiler's clamp, exactly as
`tests_ck_sweep_cost` §5 found — and no VGPR spill (19 SGPR spills, the same everywhere).
So the sweep's occupancy is set by the clamp and by the private segment, and the
private segment confirms section 1.6: the two-pass forms carry **less** private memory
than the form they replace (by 128-224 B/thread), because they hand `I_down` and `Cmx`
back and add nothing.

### 4.2 Cycles per second

300 cycles, `cpu time used` from the code's own end-of-run line (the same number
`bench/prof_0922` and `tests_ck_sweep_cost` quote), cycles/s = 300 / that.  One binary,
`bench/cksd_0922/athena.form`, for all four arms; jobs 11934363-66 on `apudev`.

| arm | `ck_sweep_cache` | `ck_sweep_form` | `cpu time used` | cycles/s | vs the default |
| --- | --- | --- | --- | --- | --- |
| four-pass, cache off | 0 | 0 | 29.701 s | 10.10 | 0.96x |
| four-pass + cache (the shipped default) | 2 | 0 | 28.445 s | 10.55 | 1.00x |
| **tm** | 2 | 1 | **20.899 s** | **14.35** | **1.361x** |
| **sd** | 2 | 2 | **22.483 s** | **13.34** | **1.265x** |

and, for scale, `tests_ck_sweep_cost` §3.2 measured the PLANE-PARALLEL kernel on the same
restart, input and nodes (a different binary, job 11933619) at 21.859 s = **13.72
cycles/s**.  So:

* **tm removes the whole cost of the spherical form.**  The 1.32x penalty
  `ck_spherical` used to carry (13.72 -> 10.10) is gone: at 14.35 cycles/s the exact
  spherical sweep is 1.046x FASTER than the plane-parallel four-pass reference, i.e. the
  geometry is now free and the restructuring pays a little on top.  (Cross-binary, so
  read that last 4.6 % as "no penalty" rather than as a speed-up.)
* **tm beats sd by 1.076x**, which is the divide count of section 1.2 showing up exactly
  where it was predicted: sd pays one divide per half layer (two per cell) to make the
  face free, tm pays one divide per face to make the layer free, and there are two half
  layers per face.  Elegance costs 7.6 %.
* The column cache is worth much less to the two-pass forms than to the four-pass one
  (there are two crossings to amortise over instead of four), so `ck_sweep_cache = 2`
  is kept as their only instantiation rather than re-swept.

**The verdict: `ck_sweep_form = 1` (tm).**  It is exact where the default is
`O(beta^2)`, 1.36x faster than the default, 1.08x faster than sd, uses 224 B/thread less
private memory than the default and the same registers, and it makes the spherical
correlated-k sweep cost what the plane-parallel one costs.

## 4.3 What remains

* `ck_implicit` is refused: the column solve's tridiagonal is the linearisation of the
  four-pass recurrence.  The Riccati form's Jacobian is a different (and simpler)
  object -- `dR/dB` and `dSc/dB` satisfy their own two-term recurrences -- but it has not
  been derived, so the implicit path still runs the four passes.
* `ck_spherical = false` is refused rather than routed through the same code.  With
  `beta = 0` the tm recurrence collapses to `R = 0`, `Sc` = the ordinary up-sweep, so the
  plane-parallel path COULD be unified with this one; it would cost a per-face `u_below`
  column and it would not be bitwise, so it was left alone.
* The 2.6 % the probe was getting wrong is a change to a production answer.  Whether the
  cs deep-hot-Jupiter production run should be restarted on the exact form is the user's
  call; nothing here changes a default.

## 5. Reproducing

```bash
cd /viper/u2/jinma/ATHENAK/athenak
cmake -B build_cksd -D PROBLEM=deep_hot_jupiter_rt && (cd build_cksd && make -j)
cd tests_ck_sweep_form
# (b) the three columns, three forms
B=../build_cksd/src/athena
for f in 0 1 2; do for c in "m0k5 0 5" "m1k3 1 3" "night 0 2"; do set -- $c
  mkdir -p col_${1}_f$f && (cd col_${1}_f$f && $B -i ../dhj_ck_sweep_form.athinput \
    time/nlim=1 problem/ck_sweep_form=$f problem/ck_dump_file=col.txt \
    problem/ck_dump_m=$2 problem/ck_dump_k=$3 \
    output1/dt=1e30 output2/dt=1e30 output3/dt=1e30 output4/dt=1e30 > run.log 2>&1)
done; done
python3 coldiff.py col_m0k5_f0/col.txt col_m0k5_f1/col.txt
# (c)
python3 ../tests_ck_sph/budget2.py col_m0k5_f1/col.txt
```
