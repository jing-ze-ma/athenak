# NOTE 2026-10-19 (Raven coordinator) — accretor: switch to accwall-1019 4d2d5628

USER approved: "go with (a) in link 2, apply the booking fix". This supersedes the accretor line in NOTE-2026-10-19-raven-prod-ab758e94.md.

- Branch `accwall-1019` @ 4d2d5628 (pushed; = ab758e94 + one commit).
- Changes:
  1. `problem/theta_wall_closed` (opt-in; needs theta_bc=hse). Every x2 flux on the theta-wall faces is zeroed except IM2 (wall pressure). This removes the hse-wall mass/energy leak: in Raven link 1 the whole Mout of -207 came through the theta walls, about 1.8x the stream inflow.
  2. Ebnd/Eirr booking fix (always on, history only): the radiation stage resets go into a plain non-RK sum `erinj_` (registered in dt_retry). The energy budget residual drops from 23 % to 1e-3 of |Ebnd|.
- Default path is BITWISE vs ab758e94 (gbit bsg/e1/hg; accretor without the key except user.hst Ebnd/Eirr).
- With the key: determinism, forced retry, restart all BITWISE on A100 and MI300A (r71).
- Deck: `acc3d_prod_v11amw.athinput` = v11am + `theta_wall_closed = true` (md5 a87998f5). Copy from Raven /raven/ptmp/jinma/acc3d_1017/prod/.
- Raven production uses this from link 2 (31291344).

ACTION (Caltech accretor copy, if running):
- Build 4d2d5628 (H200, nofma as before).
- At the NEXT link boundary, switch BIN + deck to v11amw.
- Report the md5 and the first-link Mout / mass-budget residual.
- It's a method change (not bitwise across the boundary).
Full gate report: /raven/ptmp/jinma/acc3d_1017/accwall/NOTE.md.
