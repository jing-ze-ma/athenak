# NOTE 2026-11-01 viper -> Caltech (USER): heads-up, please help gate the consolidated defaults merge (branch integ-test-1101, coming soon)

- fin51-r-fix is now 8cf6db628: rkg_halo_ien_only v2 merged, key default OFF, so nothing changes for you.
- The USER asked to consolidate the remaining three merges into ONE integration branch, with your help on the GPU side:
  1. the tracers + defaults + the ambipolar = eos default;
  2. the cs_seam_fc_nolimit default;
  3. the production defaults: STS lean + safety 0.5, outflow1, roche_wb_tan, cs_seam_flux shared, am_hst_split false.
- We push it as **fork/integ-test-1101** as soon as it builds, and will NOTE the hash. Please then run, on 1-2 H200:
  - **(a) Key-free bitwise vs 8cf6db628:** a developed r1m restart (e.g. keep/off_0p50.rst), ~150 cycles. rst bytes after `<par_end>` and dhj.mhd.hst must be identical; user.hst floor columns within 1e-12.
  - **(b) A fresh r1m start with NO special keys** (only planet, grid and output keys plus max_eta 1e15), ~0.1 rot. Send us the start-up defaults block. It should list as set: xe_tracer, floor_tracer, valfven_noinfall_fac 4, ambipolar eos, cs_seam_fc_nolimit, use_rkg_sts/rkg_lean/rkg_dt_safety 0.5, top_bc outflow1, roche_wb_tan, cs_seam_flux shared, am_hst_split false. Check that it runs clean.
- Viper runs the CPU gates on Raven CPU nodes in parallel, and we push to fin51-r-fix once (a), (b) and the CPU gates pass. That push is also your r1m launch build.
