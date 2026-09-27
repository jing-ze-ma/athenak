---
name: gpu-build-race-flag-not-compiled
description: The 22:31 GPU build of ae300379 raced with the source edit and never compiled the rt_semi_implicit READ into the pgen object; every 09-09/10 GPU ablation arm ran the semi-implicit RT source ON despite rt_semi_implicit=false in the input. Verify a flag with strings + a startup print
metadata:
  type: feedback
---

Found 2026-09-10 by object-file forensics: build_dhj_gpu's deep_hot_jupiter_rt.cpp.o (22:31:12)
had the NEW header but the OLD .cpp (the two files were written 1.2 ms apart while hipcc was
already reading). `strings -a athena | grep -c '^rt_semi_implicit$'` = 0 on every GPU binary
in bench/cs_ablate and = 1 on the CPU binary. Consequences: ctl2/nocond/nowb/norot/cache1,
det_new, det_cyc_new all ran the semi-implicit source ON (ctl == ctl2 byte-identical for that
reason); the "GPU bisect" attribution to an rt_Em FMA-contraction effect was WRONG -- the
old-vs-new GPU difference is simply semi-implicit ON vs OFF, 1-KE 3.65e-4 at cycle 1
(1e-14 of tot-E), chaos-amplified. The branch IS active on dhj (lambda*dt ~ 1e-2 at the top).

**Why:** a subagent rebuilt the GPU binary in the same second the edits landed; make saw the
.cpp as up to date. ae300379's commit message ("byte-identical on GPU with the flag false")
was never true for that binary.

**How to apply:** after ANY rebuild that adds an input flag, (1) `touch` the edited .cpp and
rebuild, (2) check `strings -a <binary> | grep -c '^<param>$'`, (3) have the pgen PRINT the
flag at startup (deep_hot_jupiter_rt.cpp prints bc_outer_maxwell that way) and read it in the
log before trusting any arm. See [[validate-the-instrument]],
[[rt-semi-implicit-changes-dhj-answer]].
