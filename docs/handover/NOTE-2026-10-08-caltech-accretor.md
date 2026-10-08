# NOTE 2026-10-08 Caltech: accretor takeover, env13 port PASS, adiabatic tests running

From: Caltech session. Answers TASK-2026-10-08-caltech-accretor-handover. Docs only.

## Port (8cecb89a, H200)
- Binary builds/bin/athena_gpu_ry_per_accretor_8cecb89a, md5 5fd715d100ffbc2ea44a6203affff15b (CUDA 12.9, sm_90;
  DeltaAI GH200 md5 dbbef84b). Input plaskett_env13.athinput md5 313d7f0b (matches).
- Smoke 4235770: 20 cycles on 1 H200 (8 MeshBlocks on one GPU, not 4 GPUs: queue), rc 0, 0 FATAL.

| cycle | time (Caltech = DeltaAI 3334844) | dt |
|---|---|---|
| 0 | 0 | 1.918640e-06 |
| 10 | 1.900746e-05 | 1.851290e-06 |
| 20 | 3.712850e-05 | 1.786984e-06 |

  All three startup lines identical (r_meas 9.61568; r_top 9.1635, max 9.6863; c_s,don 18.679, width 0.7350,
  sigma_phi 0.05444). Production speed on 1 H200: ~33 cycles/s (s1 runs below), i.e. ~1.1-1.3 h per half orbit.

## User 10-08: go beyond isothermal (steps 1 + 2), half-orbit tests on 1 H200
env13 is ideal gas gamma 5/3, but env_t_relax(_env) = 1e-4 (~50 steps) clamps T everywhere -> effectively isothermal.
- Branch **accretor-adiab-1008** (65ec7626, on top of accretor-1007): key `problem/env_t_relax_stream` (default
  env_t_relax): relaxation time of the stream/atmosphere gas outside the envelope; 0 = adiabatic, the hot ambient keeps
  env_t_relax. Gate 4235834: default keys vs 8cecb89a, 20 cycles: hst + cycle/dt lines identical, bin data identical
  (header differs only by the new key in the parameter dump). The key must be in the input file (command-line-only
  keys are FATAL in AthenaK): env13_trs.athinput = env13 + that one line.
- Runs (tlim 0.229133 = half orbit; dir /resnick/groups/carnegie_poc/jingze/accretor_1008):
  - s1_nostream 4235835: no stream, env_t_relax_env = 0 (does the envelope hold without the clamp?), running.
  - s1_env13 4235836: env13, env_t_relax_env = 0, stream/ambient relaxed as q13_s1, running.
  - s2_env13 4235882: env13, env_t_relax_env = 0 AND env_t_relax_stream = 0 (shock heating kept), queued.
- Step 3 (physical cooling: diffusion time with electron-scattering + Kramers opacity) after these, user decides.
DeltaAI has no further action.

## UPDATE 10-08 12:15 PDT: half-orbit results (t 0.2291), all 0 FATAL
| | q13_s1 (DeltaAI, all relaxed) | s1_nostream | s1_env13 (envelope adiabatic) | s2_env13 (envelope + stream adiabatic) |
|---|---|---|---|---|
| stream mass in | 8.575 | - | 8.576 | 8.562 |
| through r_meas (accreted) | 7.884 (0.919) | 1.2e-7 | 7.751 (0.904) | 3.164 (0.370) |
| out of the domain | -5e-5 | 4e-12 | 2.2e-3 | 0.73 |
| j accreted / j_K | 0.556 | - | 0.546 | 0.388 |
| min dt | 9.58e-7 | 1.78e-6 | 5.71e-7 | 4.83e-7 |
(from ryper.user.hst cols 7, 13, 8, 14; q13 numbers from numbers_q13_s1.json)
- Envelope without the T clamp is static (no stream: mass, J, tot-E unchanged to 6 digits).
- s1 = q13 within ~2 %: relaxing only the stream/atmosphere is enough to reproduce the isothermal answer.
- s2 (no cooling at all): the impact shock heats the stream and outer envelope to sqrt(P/rho) ~100-110 km/s (median, dense gas
  at r 9-10.5; c_ph 19.4) = ~30x the photospheric T; the hot layer puffs above r_meas (dense mass at r 9.6-10.5: 7.3e4 vs
  1.3e4 in s1), only 37 % crosses r_meas in the half orbit, 0.73 leaves the domain. This is the no-cooling upper bound;
  the real photosphere cools in ~1 s (docs/dev/accretor_rhd_design.md), so step 3 (cooling) decides between s1 and s2.
