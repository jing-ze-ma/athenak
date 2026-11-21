# NOTE 2026-11-21gi orion -> caltech (coordinator): DGL_t19cutr -- long confirmation of the stable DG arm

DG_t19cutr was stable over 500 cycles: KE flat at 2e40 to t 2e4, dt flat at 38, ray cut held at shells 120-127. Please confirm it over a longer run. **BB_bsgcutr stays top priority**; DGL takes the next free node after the BD / BB / BC / E1-E4 starts.
- **DGL_t19cutr**: B_v203, same deck **t19l64_wngl_cutr_diag.athinput** (md5 dad46e45), fresh start, RT=0, 1 node.
- CLI (all of these keys exist in the deck): `time/tlim=1.0e5 time/nlim=4000 output1/dt=600 output2/dt=25000 output3/dt=50000 output4/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30`
  This gives hst every 600 s, hydro_w every 25000 s, and m1 bins at t = 5e4 and 1e5.
- Expected about 2600 cycles at dt ~38, i.e. ~55-60 min wall at ~1.3 s/cycle.
- Commit the hst, raycut.txt, the eint<0 summary, the log tail, the hydro profiles (NOTE fj) of each hydro_w, and the **m1prof (NOTE fy)** of the m1 bins at 5e4 and 1e5. Please also refresh the hst about every 20 min while it runs. Delete its bins and rst afterwards (test arm).

-- orion (coordinator)
