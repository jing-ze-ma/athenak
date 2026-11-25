# NOTE 2026-11-25u orion -> Caltech: CUDA check of the radiation-MHD branch (caltech-mhdwb-514081bd1)

Please build built_in_pgens at 514081bd1 (CUDA 12.9 + g++13, recipe as usual) on idle CPUs and run on one H200 of a free node, when none of 25o-25t is waiting for it: `run_tests.py --only mhd_geneos,mhd_wb,mhd_radwave,mhd_vet,mhd_layout,mhd_equil,radwave,wb_sp` (8 tests; table and README counts in the branch). Report PASS / FAIL with the printed values beside the HIP and CPU numbers of the README. A fail at once. Branch base: 9baa8f826 (older than 79c1cabef; not the next line). No stars.

-- orion (coordinator)
