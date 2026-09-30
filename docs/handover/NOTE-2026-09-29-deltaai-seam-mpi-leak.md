# NOTE from DeltaAI: cubed-sphere seam flux MPI request leak (fixed, merged 793e03c3)

**DeltaAI, 2026-09-29 ~19:25 CDT.** User decision: merged into rt-integration. It affects every cubed-sphere run
(hydro and MHD) with more than one rank. It kills runs on Cray MPICH, and may leak memory silently under other MPIs.

## Symptom
After the c2p_track fix (NOTE-2026-09-29-deltaai-c2p-track-cuda.md), three w121 MHD arms (b3_e0/e11/e12_f13,
2 ranks) each ran ~8100 cycles (~2 min) and then aborted on both ranks:
```
Assertion failed in file ../src/include/mpir_request.h at line 497: req != NULL
MPI_Isend <- MeshBoundaryValuesCC::PackAndSendCC <- Driver::ExecuteTaskList <- Driver::Execute
```
(jobs 3272853/3272946/3272953). MPICH had run out of request handles.

## Cause
`PackAndSendFluxSeamCC` (Hydro/MHD SendFlux on the cubed sphere) posts its sends into `sendbuf[].flux_req`. The only
waits on those requests were `ClearFluxSend` / `ClearFluxRecv`, and Hydro::ClearSend/ClearRecv and MHD::ClearSend/
ClearRecv called them **only with SMR/AMR**. On a uniform cubed sphere, every stage therefore overwrote the pending
send request with a new one without waiting on it: one leaked request per off-rank seam neighbour per stage.
- viper's MPI evidently tolerates this (the 1x/10x hydro productions ran for days).
- The earlier DeltaAI hydro benchmarks ran only 1200 cycles.
- **The DeltaAI hydro production (w121prod) would also have died after a few minutes.**

## Fix: 2dfe4e95 (branch fix-seam-flux-send), merged as 793e03c3
The condition `multilevel` becomes `multilevel || use_cubed_sphere` in the four Clear functions
(hydro_tasks.cpp, mhd_tasks.cpp). These are waits only; the result is bitwise unchanged. Refinement and the cubed
sphere are mutually exclusive, so flux_req is never shared between the two exchanges.
- DeltaAI binary athena_dhj_dev_2dfe4e95ff1e (source identical to 793e03c3).
- Verified with the full b3_e0_f13 arm (job 3272979, 2 GH200): past cycle 24000 at the time of writing (3x the
  crash point), dt 4.7 -> 3.96 s as the field winds up, 14.8 ms/cycle, no errors.

## Other
- Caltech (H200, CUDA) needs both fixes, i.e. rt-integration >= 793e03c3, for any dhj MHD run. It needs the seam fix
  for any long multi-rank cubed-sphere run if its MPI has a bounded request pool.
- viper's running arms (ae767d20) are unaffected in practice, since HIP/viper MPI ran the smoke fine; rebuilding is
  not required for the scan.
- DeltaAI scan: the 12 remaining arms were resubmitted with the fixed binary on ghx4 (3273060-71); one at a time
  moves to ghx4-interactive. b3_e0_f13 is running.
