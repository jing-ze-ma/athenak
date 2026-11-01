# NOTE 2026-11-01 Caltech -> viper (USER): 2 H200 reserved for the tracer GPU tests

Reply to NOTE-...-tracer-tests-heads-up (a515b5e7).
- **GPU 1:** hpc-sm-01-05 GPU 1 is free now and reserved (WASP-121b r1m is held).
- **GPU 2:** we free a second H200 when your TASK arrives, by stopping the TOI-1518b open-top test (hpc-sm-02-03 GPU 2).
- **Capacity:** one H200 holds 3 r1m smoke arms at once (~38 GB each, the setup of our fac 2/4/10 smokes). So your 4 arms + the key-off bitwise check fit on the 2 GPUs in ~1 h.
- **Turnaround:** we build tracers-test-1101 (incremental, fmad off) and run the arms as soon as the TASK lands. Results go in one NOTE, compared against p_f4 (Mfloor, dt, E_B, bands at p > 1e-5 bar, plus night eta/Rm/x_e if you name the diagnostics).
