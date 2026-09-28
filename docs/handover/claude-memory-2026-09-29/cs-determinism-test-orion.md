---
name: cs-determinism-test-orion
description: The cs MHD restart-determinism test on orion CPUs (job 193377, /orion/ptmp/jinma/Athenak/cs_determ): protocol, the reconstructed cs production input, and where the result lands
metadata:
  type: project
---

Asked for by the viper handover ([[cs-mhd-prod-nan-rot41]]): restart twice from one rst,
200 cycles, bitwise compare. Orion has no GPU, so this is the CPU (MPI+OpenMP) version;
a clean CPU result exonerates the restart STATE and the MPI/OpenMP path, not GPU atomics.

**Run dir** /orion/ptmp/jinma/Athenak/cs_determ, binary build_dhj_mpi (5e9b1273, gcc/13 +
openmpi/4.1, OpenMP, SPR). Input cs_determ.athinput = repo deep_hot_jupiter_rt_eos.athinput
(WB polytropic + rot_potential + radiative blend + ck table already default) rewritten to the
viper cs production shape, RECONSTRUCTED from memory notes since viper's copy is not in git:
cs 128x32x32, x1max 2.0556e10, refit poly stretch (-0.068392, -2.191487, 2.464818,
-1.366698), meshblock 128x8x8 (96 blocks), grav_point_mass, bbot 3 G, max_eta 1e13, STS off,
pfloor 1e-3 / dfloor 5e-13 / tfloor 200 K, hst every cycle (dt=1 s), bin+rst every 1e4 s.

**Protocol (submit.sh, one 16x7 node on p.exclusive):** A from scratch nlim=800 (final rst
written by Driver::Finalize at exactly cycle 800) -> B1, B2 restart from that rst with
time/nlim=1000 on the SAME 16x7 layout -> B3 same restart on 28x4 (different decomposition
and reduction order, expected to differ at round-off; it calibrates what "differ" means).
The job's log ends with IDENTICAL/DIFFER lines for rst, bin and hst of B1 vs B2 and B1 vs B3.

**Trap hit while writing the input:** a regex keyed on `<meshblock>` matched the literal
`<meshblock>` inside a COMMENT of the `<mesh>` block and set the mesh nx2/nx3 to 8. Anchor
block headers at line start and check the block does not cross the next header.

**RESULT (job 193377, 2026-09-07 19:05): BIT-IDENTICAL.** A ran 800 cycles (t=1.88e4 s,
3.44e6 zone-cycles/s on one node); B1, B2 (16x7) and even B3 (28x4, different decomposition
and reduction order) restarted 200 cycles to t=2.22e4 s and their final rst, bin and hst
files have the SAME md5 (rst 1a74c0c5..., bin 47248e1a..., hst 841dc5bf...). So on CPU the
restart state is exact and the cs seam exchange / EOS resistivity / WB cache / radiative
blend path is deterministic, even across rank layouts. The 5e-6 divergence seen on viper
([[cs-mhd-prod-nan-rot41]]) is therefore GPU-specific (atomics / reduction order on the
MI300A) or a difference between the two viper runs' binaries -- not the restart machinery.
Next on this thread: repeat B1/B2 on viper's GPU after 09-12 with the SAME binary.
