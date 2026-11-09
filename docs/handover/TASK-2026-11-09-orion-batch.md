# TASK 2026-11-09 viper -> Orion: CPU batch (USER: "anything that can be run on cpu can go to orion")

Run these self-contained TASKs in this order (all CPU, MPI build lines inside each). Answer each with its own NOTE
as the TASK says. Independent of each other; run in parallel if Orion has the cores.

1. `docs/handover/TASK-2026-11-09-orion-hgpos.md` - HG positivity: time2_adm_eflr + implicit_pos_transport_t gates
   (fork/hgpos-1109b b798f1a30 vs prod-1112-cand f34ada236). NOTE: NOTE-2026-11-09-orion-hgpos.md.
2. `docs/handover/TASK-2026-11-09-orion-r136fast.md` - R136 floor-gas cap + re-block identity gates
   (fork/fastfloor-1109 138e6ecb5). NOTE: NOTE-2026-11-09-orion-viper-r136fast-gates.md.
3. `docs/handover/TASK-2026-11-09-orion-raycheap.md` - physical-clock ray-average gates + 2nd-order convergence check
   (fork/raycheap2-1109 a39c9ec49 vs prod-1112-cand 3330bc0dd). NOTE: NOTE-2026-11-xx-orion-raycheap.md.

Not included: TASK-2026-11-08-orion-hg-mesa-lin.md (USER declined on Orion, see NOTE-2026-11-08-orion-viper-hg-mesa-declined.md).
4. `docs/handover/TASK-2026-11-09-orion-lintol.md` - R136 small (831x8x8) Picard stop-test accuracy: find the cheapest
   stop setting reaching ~1.5e-8 vs the reference (deck/data in docs/handover/lintol-1109/; fork/prod-1112-cand).
   NOTE as the TASK says.
5. `docs/handover/TASK-2026-11-09-orion-raycheap2.md` - restart fix gates (fork/rstabe-1109 4e9f602c8: hesdirk2
   auto_be back-off in the rst) R1-R4 + convergence reconciliation C1-C3 (24 ranks, tight tolerance). NOTE as the TASK says.
6. `docs/handover/TASK-2026-11-09-orion-hgblend.md` - HG transverse blend implicit upwind (fork/hgpos-1109b
   1f314023a): keys-absent / 1 vs 4 ranks / restart gates + active-path check. NOTE as the TASK says.
7. `docs/handover/TASK-2026-11-09-orion-raycheap3.md` - why the base VET/M1 time scheme is 1st order: time-order
   isolation arms B0 C1 C2 C3 G1 on the forced test star (fork/raycheap2-1109 9b9951606, 24 ranks). NOTE as the TASK says.

9. `docs/handover/TASK-2026-11-09-orion-raycheap4.md` - CORRECTED 11-09 (transient dumps): time-order test of the base
   VET/M1 scheme on a smoothly ramped filling transient, runs ending at t = 2.5 / 5 / 15, arms P0-P3 x 4 levels
   (fork/raycheap2-1109 ff2ed77f1, REBUILD; 24 ranks). If you already ran the first version (pulse/ decks, t_end 30),
   discard it and rerun from this text. NOTE as the TASK says.
8. `docs/handover/TASK-2026-11-09-orion-hglagt.md` - HG all-lagged transverse upwind (fork/hgpos-1109b
   7fb7af05a): keys-absent / 1 vs 4 ranks / restart gates + active path. NOTE as the TASK says.

10. `docs/handover/TASK-2026-11-09-orion-raycheap5.md` - which VET update restores 2nd order in time: arms U1-U7
   (rebuild / x2_fields 15 / lag / extrap / predict / lat offdiag none / every 2) on the pulse2 test
   (fork/raycheap2-1109 fa7103f33 = the ff2ed77f1 build, 24 ranks). NOTE as the TASK says.

11. `docs/handover/TASK-2026-11-09-orion-raycheap6.md` - gas-driven closure change: lag / extrap2 / predict / rebuild,
   with and without the tangential tensor (gasmod decks), + pulse2 re-test (fork/raycheap2-1109 f872899da, REBUILD; + (c) vet_source = gas arms;
   24 ranks). NOTE as the TASK says.

12. `docs/handover/TASK-2026-11-09-orion-raycheap7.md` - VET update order with gas_feedback on (source relaxation
   active): X9 / PR / RB / PRg / X9t / PRt, 4 levels (fork/raycheap2-1109 8532af080 = the raycheap6 build). NOTE as the TASK says.

13. `docs/handover/TASK-2026-11-09-orion-r136cap13.md` - gate of the COMBINATION prod-1113-cand bc35bfca8
   (= prod-1112-cand 4e9f602c8 + fastfloor-1109 floor cap): keys absent vs 4e9f602c8 / cap 1 vs 4 ranks / cap restart
   (incl. auto_be back-off) / cap counters fire. NOTE as the TASK says.

14. `docs/handover/TASK-2026-11-09-orion-raycheap8.md` - evolving-gas time-order test of the VET update (equilibrium IC,
   gas feedback, stability pre-check inside) + CPU check that the stage-2 rebuild reaches the operator. NOTE as the TASK says.

15. `docs/handover/TASK-2026-11-09-orion-hgfct.md` - HG FCT limiter on the lagged transverse remainder (fork/hgpos-1109b
   a4bac617e): keys-absent vs 7fb7af05a / 1 vs 4 ranks / restart + forced active path (alpha 1e-3). NOTE as the TASK says.

16. `docs/handover/TASK-2026-11-09-orion-raycheap9.md` - isolate the common O(dt) term on Ar4: hydro off / Strang /
   alternate / radiation force off / no clump / be, base PR, 4 levels (reuse the raycheap8 build and Ar4 decks). NOTE as the TASK says.

17. `docs/handover/TASK-2026-11-09-orion-raycheap10.md` - R alone on hesdirk2 (new dbg_hydro_skip_t2, fork/raycheap2-1109
   432b5055c, REBUILD) on Ar4, then R-internal toggles if still 1st order. NOTE as the TASK says.

18. `docs/handover/TASK-2026-11-09-orion-raycheap11.md` - R alone 1st order: stage order of the exchange (frozen closure,
   no-VET, Newton opacity, coupling off, TR-BDF2 arms) on Ar4 (fork/raycheap2-1109 addb642d0, REBUILD). NOTE as the TASK says.

19. `docs/handover/TASK-2026-11-09-orion-raycheap12.md` - time order of the M1 scheme on clean 1-D benchmarks (EQUIL
   non-stiff/stiff, diffusion pulse, Marshak) x hesdirk2 / trbdf2 / be x m1|eddington / vet_col. NOTE as the TASK says.

Already answered (skip): 1 hgpos, 2 r136fast, 3 raycheap.
