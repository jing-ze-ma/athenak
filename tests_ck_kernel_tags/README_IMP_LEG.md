# IMP and LEG: ck_implicit and rt_layer_legacy as compile-time tags on rt_chain_ck (2026-09-22)

Follow-up to section 5 of `README.md`.  **Verdict: correct and bitwise on CPU, but no
measurable gain for the production default.**  The production sweep (`ck_spherical`,
`ck_sweep_form = tm`) runs at the same speed as before.  The plane-parallel sweep runs
15 % faster.  The change is not committed.  The patch is
`bench/tags_imp_leg_0922/imp_leg.patch` (against HEAD `3b6e0e52`, 1 file,
`src/utils/two_stream_rt.hpp`).

## 1. The change

`launch_ck_chain` gets two more tags next to `FOP`: `IMP` (`ck_implicit`) and `LEG`
(`rt_layer_legacy`).  Inside the kernel:
`if (IMP && ckskip_ && ...) return;`, `const bool jck = IMP && ckjacp_;`, `if (LEG)`
in place of `if (layer_legacy)`.  They are dispatched in `launch_ck_cache` and nested so
that only the combinations the startup guards allow get instantiated:

| branch | FOP | IMP | LEG | per tier |
| --- | --- | --- | --- | --- |
| `ck_implicit` (+ `ck_impl_frozen_op`) | 0/1 | 1 | 0 | 24 |
| `rt_layer_legacy` (refused with SPH, BSP and ck_implicit, so only SPH = BSP = 0) | 0 | 0 | 1 | 3 |
| **production / everything else** | 0 | 0 | 0 | 12 |
| `ck_sweep_form` tm/sd (refused under ck_implicit; legacy is refused with SPH) | 0 | 0 | 0 | 4 |

The kernel count goes from 28 to 43 per tier (112 to 172 in total).  The two HIP builds
ran side by side and finished 17 s apart, so the extra compile time is small.
cpplint: no new note (the pre-existing `whitespace/blank_line` note is still there).

## 2. Code object (tier 136, gfx942)

These counts cover the kernel plus every device function it calls
(`bench/tags_imp_leg_0922/count2.py`).  The old `count.py` counted only the outlined
`operator()` bodies, and in the patched build most chain bodies are inlined into the
kernel, so the two methods cannot be compared.  Both builds were counted with the new
method.  VGPR use is 128 in every instantiation, before and after.

| instantiation | HEAD instr / scratch ops | IMP=LEG=0 | change |
| --- | --- | --- | --- |
| SPH1 BSP0 CCH2 **FRM1 (production, tm)** | 25470 / 1580 | **18840 / 1329** | -26 % / -16 % |
| SPH1 BSP0 CCH2 FRM0 (four-pass) | 34424 / 2390 | 22437 / 1332 | -35 % / -44 % |
| SPH0 BSP0 CCH1 (plane-parallel default) | 27581 / 1922 | 17665 / 1208 | -36 % / -37 % |

Caveat: the metadata reports `vgpr_spill_count` 0 for the HEAD production kernel and
407 for the patched one.  This is probably an accounting artefact: HEAD's spills sit in
its outlined callee, and the kernel-level metadata does not count them.  The private
segment is 27104 bytes in HEAD and 27024 in the patch.

## 3. GPU timing: apudev job 11942108 (8 arms, 3:47 wall)

The recipe is the same as in README section 4: restart `cs_mhd_prod3/rst/dhj.00567.rst`,
the production input, 2 ranks / 2 GPUs, 300 cycles.  HEAD and new arms alternate.

| setting | HEAD `cpu time used` | new | cycles/s HEAD to new |
| --- | --- | --- | --- |
| defaults (spherical, tm) | 20.28 / 19.88 / 20.44 (mean 20.20) | 20.46 / 19.97 / 20.05 (mean 20.16) | 14.85 to 14.88 (**+0.2 %, noise**) |
| `ck_spherical = false` | 21.02 | 17.80 | 14.27 to 16.86 (**-15.3 % time**, 1 sample) |

The production tm kernel is 26 % smaller but not faster.  As in README section 3, the
instruction count does not predict the time here.

## 4. Bitwise

**CPU** (`bench/tags_imp_leg_0922/gate/gate.sh`, `gate_leg.sh`; 20 cycles; the `.hst`,
22 bin dumps and 2 restarts per arm): every arm is **BITWISE**.  The arms were: defaults
(sph off), `ck_spherical`, `ck_implicit`, `ck_implicit` + `frozen_op`, `ck_implicit` +
`frozen_op` + `ck_spherical`, `ck_implicit` + `reuse_jac = 1`, and `rt_layer_legacy`.
`rt_layer_legacy` cannot be set on the command line (the parameter does not exist in the
input), so that arm uses `gate/dhj_leg.athinput`.  Its `.hst` differs from the defaults
arm, which shows the switch is live.  Checksums are in `gate/md5.txt`.

**GPU** (the timing runs): the `.hst` is identical.  The bin dump and the restart differ
at float32 round-off (at most 4.8e-6 relative, in about 200 of 786k cells per variable;
`dens` is identical).  Each binary is deterministic on its own (h1 = h2 = h3,
n1 = n2 = n3).  The FOP commit behaves the same way: its r0 and r1 outputs in
`bench/cktag_0922/gpu` also differ.  So "bitwise" in README.md holds for CPU only; on
GPU, a change in code generation moves the last bits.
