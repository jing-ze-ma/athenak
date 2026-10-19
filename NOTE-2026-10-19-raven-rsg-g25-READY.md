# NOTE 2026-10-19 raven -> caltech: RSG g25 production READY (reply to NOTE 8c76c3c7, TASK bd1e651d)

**READY.** The Raven reduced smoke is clean on the gated deck. Raven chain 31285583 -> 84 -> 85 (6 nodes x 4 A100) has READY=1.

## Gate
- Code: `origin/fin51-r-fix` tip **a2713be0** (unchanged at gate time). Raven A100 binary md5 4d23e8c9e709ee60e67a035dcb362ceb.
  Viper MI300A rocm/7.1 build of the same tree: md5 f8c8e8650696769c1e4467dfb89803a0.
- Deck `rsg-g25-prod/rsg_g25_prod.athinput` md5 **b49f142df8fc4a4a73694989dbf634cd** (IC pass3b = seam ON, fixed point R_ph 5.80481e13 cm, 2 passes, last dR 0.034 cells).
- Smoke (cs N=48, 24 blocks 591x24x24, 200 cycles to t 6.57e4 s), on Raven 4 A100 and on Viper 2 MI300A:
  - C and R: rc 0, no FATAL/NaN/DIVERGED/DT_RETRY, Picard NON-CONVERGED=0 (mean 2.96, max 5).
  - Restart R vs C: BITWISE (hst rows and bin data).
  - Run-to-run C vs C2: BITWISE (hst and all bins).
  - Peak memory 29.5 GB per A100 at 2.0M cells per GPU.

## Files in `docs/handover/rsg-g25-prod/` (md5)
- rsg_g25_prod.athinput b49f142df8fc4a4a73694989dbf634cd: Raven paths and the Raven layout, 24 MeshBlocks 591x64x64.
- ic_rsg_g25.txt 87918742791146d42167cd43eaa66e04: `problem/he_ic_file`.
- rosseland_mesakap_4c0e8b64d26a.txt fa0f4f420022bd6c0fec4eb1803619b0: `problem/he_opac_table`.
- planck_ext2_x0.653848_z0.019757_d38e5526afef.txt a9a8d6d31bed69322fe038bba2a40377: `problem/he_planck_table`.
- diag_rsg_g25.txt, prov_rsg_g25.json: tool diagnostics and provenance, for reference only.
- run.cfg values (Raven): `TLIM=1.0e9`, `XKEYS=""`, `STOP_NC=1`. Set BIN and BIN_MD5 to your a2713be0 H200 build.

## Caltech: what to do
1. Fix the three paths (he_ic_file, he_opac_table, he_planck_table) to your copies. Change nothing else.
   - The 24 x 591x64x64 MeshBlocks stay as they are. Your layout is 3 blocks per GPU on 8 H200, set by the rank count.
2. Release links **4993305 -> 4993306 -> 4993307** with READY=1 only now, on this READY.
   - If Raven later HOLDS 31285583, hold yours too, and keep them held until a new NOTE says otherwise.
3. Race rule: the first site to START wins. The loser is cancelled only after the winner has started.
4. **Never mix sites within one run.** The layout differs (Raven 1 block per A100, Caltech 3 per H200), so a restart across
   sites is not bitwise. A run continues only from its own site's rst.
5. Restart equivalence was checked only on the reduced layout. Check the first production link once with a restart.
6. Report: the md5 of your build, the start time, s/cycle, memory per GPU and the Picard summary.

## Background check (smoke C, t = 0 -> 6.57e4 s; R_ph 5.80481e13)
- 1.0-1.08 R_ph: max |dln rho| 1.3e-3.
- 1.1-3 R_ph (floor):
  - |dln rho| <= 1.1e-2.
  - Infall only at the sponge terminal speed (g/rate ~ 0.06 km/s; at most 0.2 km/s just above 1.1 R_ph). No cell is below -1 km/s.
  - The sponge is active on all 110 floor cells (rho < 2e-17).
- One startup artefact at the floor junction (1.085-1.11 R_ph):
  - The atm-force module (phi 68 at the kappa jump) followed by the equil re-march leaves a 10x density bump (1.5e-16 over the
    floor 1e-17) above a one-cell trough.
  - That cell refills at -3 km/s; the change is dln rho 1.04 in that one cell only.
  - The mass involved is ~1e18-1e19 g (rho is 7 decades below the photosphere). Not a blocker; watch it in the first days.
