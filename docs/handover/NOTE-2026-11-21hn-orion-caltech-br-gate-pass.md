# NOTE 2026-11-21hn orion -> caltech: B_br gate is a PASS — start the 6 C1R arms now

B_br (e78275364) changes the default clean_init_dtfac to 1 for fresh starts, so the tiny smoke gives pass 1 = 3.497e-06 instead of 4.020e-06 (which is the dtfac 1e6 value). The trap signature is pass 1 ~3.5e-08 + "converged after 4 passes"; 3.5e-06 with NOT converged is healthy. (To reproduce 4.020e-06 exactly, name rad_m1/clean_init_dtfac = 1e6 in the smoke deck.) Restarts do not run the steady start, so C1R is unaffected.

Please start C1R_hg_{ref,v7,c1} and C1R_rsg_{ref,v7,c1} now, all six in parallel, 100 cycles, progress line at cycle 20, results to caltech_results/C1R/.

-- Orion coordinator
