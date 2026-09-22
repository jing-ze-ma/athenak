# T0 brief (launch with subagent_type "worker"; user order 2026-09-22)

Run phase T0 (go/no-go gate, NO code changes) from tests_ck_implicit/DESIGN_tm.md section 4
(~lines 141-160; read sections 0 and 4 first). Deliverable: tests_ck_implicit/T0_RESULTS.md with
GO or NO-GO and the numbers. Production configuration: ck_spherical = true, ck_beam_sph = true,
ck_sweep_form = 1 (tm), semi-implicit gas coupling.

A. Existing dumps, no job: T(p) difference above 0.01 bar (day side, night side, global;
   1e-6..1e-2 bar) between bench/tm_prof_growth/g_cfl03 and g_cfl015 at equal simulation time;
   max and rms |dT| per pressure level.
B. One apudev job (2 GPUs, 1 node, < 15 min) in bench/impl_t0_0922/: restart
   bench/cs_mhd_prod3/rst/dhj.00567.rst with the prod4 <problem> block as
   bench/tm_prof_growth/prof_tm did, binary bench/cs_mhd_prod4/athena (touch nothing else in that
   directory; prod4 jobs 11941995-8 are the production). A command-line override only works for
   keys present in the input file: put new keys in the input. Turn on the per-cell RT report
   (rt_cell_report / rt_report_every, src/pgen/deep_hot_jupiter_rt.cpp:628); record the
   sweep-to-gas gap (rt_desum) per MeshBlock and the rt_de_max clip count past the restart
   transient.
Criterion: GO only if |gap| > ~1 % of the column source at production dt, or upper-atmosphere T
moves > ~5 K between CFL 0.3 and 0.15. NO-GO otherwise. If a measurement cannot be made as
specified, say so and state what was measured instead. Final reply under 200 words.
Note: the viper file (inode) quota was full on 09-22; check headroom before big snapshots.
