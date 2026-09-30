---
name: red-giant-common-dt-collapse-571e5
description: "COMMON DEATH at 5.71e5 (09-10 18:00): V9f, V10 AND old-RT V9_dense all die of the EXPLICIT radiative-conduction runaway in an underdense supersonic downdraft column at 1.058 R on a cube EDGE, 20-30 cells BELOW the photosphere; nothing to do with the corona/join/top; no scheduled event; the conduction dt is evaluated from the previous state so a 4.9 s step is taken when the limit is already 1e-4 s"
metadata:
  type: project
---
Supersedes the corona/open-top suspicion for the 5.7e5 deaths. Verified from the
`### dt COLLAPSE` diagnostic in out.txt (mesh.cpp:970-1005): limiter = conduction in all three
runs, hydro dt stayed 30.6 s.

**Where**: V9f (0,2,2,315) r=3.3876e12 = 1.058 R, block 90, panel edge 7 cells from the vertex;
V10 (0,2,2,313) same column; V9_dense block 69 panel-4 edge. Photosphere is 3.4464e12 (1.077 R),
the corona join 3.65e12 -> the kill is INSIDE the star; V9f/V10/V9_dense differ only above 3.65e12
and their hst tot-E agrees to 5 digits, i.e. the surface convection is identical in all three.

**What**: a supersonic convective downdraft (v_r -1.24e6, Mach 2.5) evacuates its column: rho falls
11x between 5.2e5 and 5.6e5 to 5.7e-10 (62x below the shell median), T 2752 K vs 9625 K median.
Census at i=315: n(v_r<-1e6) 0 -> 219 -> 742 -> 1138 -> 1425 (t=4.6e5..5.6e5), grid-locked
4-down/4-up pattern at the 24 vertex columns, identical to 2-3 digits in V9f and V10.
Then explicit radiative conduction into the cell with 1/100 the heat capacity of its neighbours:
T 2752 -> 26145 -> 8.3e6 -> 3.16e10 (table ceiling) over ~500 cycles then ONE cycle; kappa_rad
1.5e11 -> 4.9e28; dt_cond ~ T^-3. dt is computed from the PREVIOUS state, so a 4.87 s step
is taken when the limit is already 6.6e-5 s (V9f c19783->19784). tot-E x35 (V9f), x5600 (V10).
Flux limiter saturated (s1=55). NOT cold collapse, NOT dfloor, NOT the RT rescue messages
(those are at i=385-392 in the atmosphere and start AFTER the dt decline began at c18500).

**Why 5.71e5**: no scheduled event (hst 1e4, bin 2e4, rst 2.5e5, mlt_alpha=0, ck_dump_t2=1e6,
sponge off); it is the growth time of the surface downdraft census.

**Next test (agent's recommendation, not yet run)**: restart V9f from rst at 5.0e5 with the
radial conduction sub-cycled/implicit at its own local limit (or re-evaluate dt_cond AFTER the
operator / cap kappa_rad by free-streaming before NewTimeStep); control: cfl_number=0.1 should
only postpone. If it passes 5.72e5 the explicit operator is the killer and the underdense
downdraft is a resolution problem. Related: [[red-giant-explicit-conduction-explosion]],
[[red-giant-vertex-chimney]], [[red-giant-dt-collapse-solved]] (the earlier 1.9e5 collapse was
the same operator at the photosphere, cured then by the tau blend).
Scripts: /orion/ptmp/jinma/Athenak/red_giant/_analysis_0910/task_dt{collapse,2,3,4,5,6}.py.
