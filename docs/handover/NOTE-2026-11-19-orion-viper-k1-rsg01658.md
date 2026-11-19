# NOTE 2026-11-19 orion -> viper (USER): staggered K1 on rsg.01658 holds dt 1000 from the restart (your Q4)

Result:
- Run: viper pHG262 13199032, 4 nodes = 8 MI300A, Orion step coord_k1_rsg01658, 2026-11-19 09:03 CET.
- Binary: Orion K line 1645edc41 (klwb2, md5 53d5ae254d3b8cc2097fecbeb87aabff). This is the staggered K1 (F solve) with the
  WB polytropic hybrid.
- Keys: our k1_rsg.args (prod K1 keys; implicit_fb_local and the KE cap are not in this line).
- Restart: rsg.01658 (your file, stored cycle 750329), nlim 750333, no restart output.

| cycle | dt | elapsed (s) |
|---|---|---|
| 750329 | 1000 | 128.7 (includes start-up) |
| 750330 | 1000 | 147.9 |
| 750331 | 1000 | 157.4 |
| 750332 | 1000 | 163.7 |
| 750333 | 1000 | 176.5 |

- No dt cut, no dt_retry, no backward-Euler step, no stage fallback.
- 0 CLEAN WARNING, 0 FATAL.
- Floor clips 0: neither "solved E <= e_floor" nor "written-back eint <= 0" fired in any solve.
- Gas Newton fallbacks 0.
- Steady speed ~9.6 s/cycle on 8 MI300A (this binary predates our restart-I/O and memory work).

So on the same restart where your clean-cc cell layout needs dt 250 after two cuts, the staggered K1 takes dt 1000 from its
first step. Relevant differences (your Q1/Q4):
- The restart imports prod's hesdirk2 K1 slope (t2 slope import).
- The F solve has the drag b = beta/(1+beta).
- The face flux is upwinded.
- The exchange is the K-line form.

Our first-iterate exchange is the same inconsistent form as yours, unless clean_xt_exact=cons. This run did NOT use cons, so the
first-iterate issue alone does not explain your failure.

Run dir on viper: /viper/ptmp2/jinma/orion_stag2/vgpu/runs/coord_k1_rsg01658 (cmd.txt, run.log).

-- orion
