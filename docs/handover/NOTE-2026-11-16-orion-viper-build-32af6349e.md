# NOTE 2026-11-16 Orion -> viper (USER): build origin/clean-cc 32af6349e for the GPU tests (start G5 now)

origin/clean-cc = **32af6349e** (push 10). It adds to 51169613a / clean-cc-int 4680ce500 only the
GPU-parallel FV ray sweep: cells are grouped into upwind dependency levels per direction, and each sweep
pass is one launch per level over all cells x directions (parallel width ~28x larger on an sp wedge,
~350x on a Cartesian 32^3 block). It is bitwise identical to the sequential sweep on every regression
output (CPU), so all references made on 4680ce500 stay valid.

Please build 32af6349e (build_cuda_raven.sh from TASK-2026-11-16-viper-cc_gpu1-a100) and start the G5
timing now; it needs no references. The final TASK files for stage 1 and for the decisions package (both
naming 32af6349e) follow within the hour; G1-G4 and the decision arms can wait for them.
