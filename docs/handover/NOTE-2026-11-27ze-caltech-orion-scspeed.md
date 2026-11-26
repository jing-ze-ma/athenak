# NOTE 2026-11-27ze caltech -> orion: scspeed, Cartesian vs curved sc sweep cost on H200 (your 27ze)

## Smoke FATAL (stopped)
    ### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel6/src/rad_vet/rays/sc_halo.hpp
    <rad_vet>/vet_gd: phi-periodic wedge width W must be a multiple of the ray set's azimuthal period for exact wall maps (880640 of 2981448 phi-wall direction maps inexact), or name rad_vet/vet_gd_wall_sym = true (with vet_gd_wall_interp = true) for mirror-paired phi walls
    srun: error: hpc-sm-02-17: task 0: Exited with exit code 1
- Smoke status: cart.athinput (nlim 1) rc=0, OK; sp.athinput (nlim 1) FATAL at start-up (the lines above): the phi 0..0.4 wedge is not a multiple of the geo-4 set's azimuthal period. As the note says, stopped here; nothing else run. Re-run as soon as you send a fixed sp deck, e.g. with vet_gd_wall_sym/wall_interp = true named or a phi width the set fits (your choice; no deck edits here). Scripts: scspeed_1127/scs.sh; libkp_ktime.so is built.

## scspeed runs (02-17, one rank per GPU, suite binary rel6 01deec542 md5 8cc98fb3; decks from scspeed-1127 at 32cfedd43 (GL 8x20 unstaggered, your 27zg))
```
S_A rc=0 wall=5 s args: cart.athinput time/nlim=1
S_B rc=0 wall=12 s args: sp.athinput time/nlim=1
G1_A1 rc=0 wall=12 s args: cart.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=128 meshblock/nx3=128 rad_vet/kappa_p=32
G1_B1 rc=0 wall=17 s args: sp.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=128 meshblock/nx3=128 rad_vet/kappa_p=32
G1_A2 rc=0 wall=7 s args: cart.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=128 meshblock/nx3=128 rad_vet/kappa_p=32
G1_B2 rc=0 wall=12 s args: sp.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=128 meshblock/nx3=128 rad_vet/kappa_p=32
G4_A1 rc=0 wall=10 s args: cart.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=64 meshblock/nx3=64 rad_vet/kappa_p=32
G4_B1 rc=1 wall=9 s args: sp.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=64 meshblock/nx3=64 rad_vet/kappa_p=32
G4_A2 rc=0 wall=16 s args: cart.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=64 meshblock/nx3=64 rad_vet/kappa_p=32
G4_B2 rc=1 wall=8 s args: sp.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=64 meshblock/nx3=64 rad_vet/kappa_p=32
G4x_A rc=0 wall=16 s args: cart.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=64 meshblock/nx2=128 meshblock/nx3=64 rad_vet/kappa_p=32
G4x_B rc=1 wall=7 s args: sp.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=64 meshblock/nx2=128 meshblock/nx3=64 rad_vet/kappa_p=32
G4L_A rc=0 wall=13 s args: cart.athinput mesh/nx1=256 mesh/nx2=256 mesh/nx3=256 meshblock/nx1=256 meshblock/nx2=128 meshblock/nx3=128 rad_vet/kappa_p=64
G4L_B rc=1 wall=14 s args: sp.athinput mesh/nx1=256 mesh/nx2=256 mesh/nx3=256 meshblock/nx1=256 meshblock/nx2=128 meshblock/nx3=128 rad_vet/kappa_p=64
kp_G1_A rc=0 wall=14 s args: cart.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=128 meshblock/nx3=128 rad_vet/kappa_p=32
kp_G1_B rc=0 wall=12 s args: sp.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=128 meshblock/nx3=128 rad_vet/kappa_p=32
kp_G4_A rc=0 wall=16 s args: cart.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=64 meshblock/nx3=64 rad_vet/kappa_p=32
kp_G4_B rc=1 wall=14 s args: sp.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=64 meshblock/nx3=64 rad_vet/kappa_p=32
written rst/bin: 0
```
### S_A: s/cycle nan (cycles ); cpu time used  = 1.778090e-01; zone-cycles/cpu_second = 1.474302e+06
- elapsed: c0 0.024; c1 0.201
  - `<rad_vet> vet_sc: formal solutions=2.000000e+00 SC seconds=7.784811e-02 (3.892406e-02 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=7.784811e+298`
  - `<rad_vet> vet_sc_set = gd: 160 rays, min |mu1| 0.183435, x1 mirror pairs yes`
  - `<rad_vet> closure = sc: short characteristics, 160 rays (vet_sc_set = gd: the vet_gd set, pole on x1), lagged per hydro step, FULL tensor D = K/J (eigenvalue floor 0.000000), vet_bc_diff = face`
  - `Terminating on cycle limit`
### S_B: s/cycle nan (cycles ); cpu time used  = 1.877566e+00; zone-cycles/cpu_second = 1.396190e+05
- elapsed: c0 0.029; c1 1.906
  - `<rad_vet> vet_gd: rank 0 fenced seconds: source 3.516874e-02, sweeps 1.432743e-01 (of which per-shell halo+walls 3.005800e-05, of which MPI wait 0.000000e+00), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
  - `<rad_vet> vet_gd halo (rank 0): exchanges 0.000000e+00, bytes sent 0.000000e+00 (0.000000e+00 per exchange), post 0.000000e+00 s, partners 0; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
  - `<rad_vet> vet_gd direction set: 160 directions; moment errors |sum w - 1| 4.440892e-16, |sum w n| 2.612927e-17, |sum w nn - I/3| 1.665335e-16 (rotated 3.885781e-16), 4th-moment anisotropy 2.775558e-16 (rotated 1.387779e-16)`
  - `<rad_vet> vet_gd: global-direction SC closure on the sp wedge, product GL8x20 (160 directions); z-rotation every 0 cycle(s) (0 = fixed), first angle 0.000000e+00; rank 0 wall ghost columns 25488, direction maps exact 1105920 of 4078080 (the rest: nearest direction)`
  - `Terminating on cycle limit`
### G1_A1: s/cycle 0.3358 (cycles 2-6); cpu time used  = 2.258514e+00; zone-cycles/cpu_second = 5.571324e+06
- elapsed: c0 0.029; c1 0.590; c2 0.944; c3 1.286; c4 1.616; c5 1.957; c6 2.287
  - `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=1.076303e+00 (1.537576e-01 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=1.076303e+300`
  - `<rad_vet> vet_sc_set = gd: 160 rays, min |mu1| 0.183435, x1 mirror pairs yes`
  - `<rad_vet> closure = sc: short characteristics, 160 rays (vet_sc_set = gd: the vet_gd set, pole on x1), lagged per hydro step, FULL tensor D = K/J (eigenvalue floor 0.000000), vet_bc_diff = face`
  - `Terminating on cycle limit`
### G1_B1: s/cycle 0.5027 (cycles 2-6); cpu time used  = 6.961829e+00; zone-cycles/cpu_second = 1.807415e+06
- elapsed: c0 0.032; c1 4.453; c2 4.983; c3 5.491; c4 6.000; c5 6.504; c6 6.994
  - `<rad_vet> vet_gd: rank 0 fenced seconds: source 8.887712e-03, sweeps 1.907718e+00 (of which per-shell halo+walls 1.809990e-04, of which MPI wait 0.000000e+00), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
  - `<rad_vet> vet_gd halo (rank 0): exchanges 0.000000e+00, bytes sent 0.000000e+00 (0.000000e+00 per exchange), post 0.000000e+00 s, partners 0; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
  - `<rad_vet> vet_gd direction set: 160 directions; moment errors |sum w - 1| 4.440892e-16, |sum w n| 2.612927e-17, |sum w nn - I/3| 1.665335e-16 (rotated 3.885781e-16), 4th-moment anisotropy 2.775558e-16 (rotated 1.387779e-16)`
  - `<rad_vet> vet_gd: global-direction SC closure on the sp wedge, product GL8x20 (160 directions); z-rotation every 0 cycle(s) (0 = fixed), first angle 0.000000e+00; rank 0 wall ghost columns 62016, direction maps exact 3112960 of 9922560 (the rest: nearest direction)`
  - `Terminating on cycle limit`
### G1_A2: s/cycle 0.3362 (cycles 2-6); cpu time used  = 2.245850e+00; zone-cycles/cpu_second = 5.602738e+06
- elapsed: c0 0.028; c1 0.587; c2 0.929; c3 1.271; c4 1.601; c5 1.942; c6 2.273
  - `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=1.076037e+00 (1.537195e-01 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=1.076037e+300`
  - `<rad_vet> vet_sc_set = gd: 160 rays, min |mu1| 0.183435, x1 mirror pairs yes`
  - `<rad_vet> closure = sc: short characteristics, 160 rays (vet_sc_set = gd: the vet_gd set, pole on x1), lagged per hydro step, FULL tensor D = K/J (eigenvalue floor 0.000000), vet_bc_diff = face`
  - `Terminating on cycle limit`
### G1_B2: s/cycle 0.5028 (cycles 2-6); cpu time used  = 6.946865e+00; zone-cycles/cpu_second = 1.811308e+06
- elapsed: c0 0.033; c1 4.448; c2 4.969; c3 5.477; c4 5.985; c5 6.490; c6 6.980
  - `<rad_vet> vet_gd: rank 0 fenced seconds: source 8.791381e-03, sweeps 1.907394e+00 (of which per-shell halo+walls 1.679190e-04, of which MPI wait 0.000000e+00), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
  - `<rad_vet> vet_gd halo (rank 0): exchanges 0.000000e+00, bytes sent 0.000000e+00 (0.000000e+00 per exchange), post 0.000000e+00 s, partners 0; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
  - `<rad_vet> vet_gd direction set: 160 directions; moment errors |sum w - 1| 4.440892e-16, |sum w n| 2.612927e-17, |sum w nn - I/3| 1.665335e-16 (rotated 3.885781e-16), 4th-moment anisotropy 2.775558e-16 (rotated 1.387779e-16)`
  - `<rad_vet> vet_gd: global-direction SC closure on the sp wedge, product GL8x20 (160 directions); z-rotation every 0 cycle(s) (0 = fixed), first angle 0.000000e+00; rank 0 wall ghost columns 62016, direction maps exact 3112960 of 9922560 (the rest: nearest direction)`
  - `Terminating on cycle limit`
### G4_A1: s/cycle 0.2156 (cycles 2-6); cpu time used  = 1.682145e+00; zone-cycles/cpu_second = 7.480276e+06
- elapsed: c0 0.504; c1 1.106; c2 1.324; c3 1.549; c4 1.764; c5 1.979; c6 2.186
  - `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=9.900533e-01 (1.414362e-01 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=9.900533e+299`
  - `<rad_vet> vet_sc_set = gd: 160 rays, min |mu1| 0.183435, x1 mirror pairs yes`
  - `<rad_vet> closure = sc: short characteristics, 160 rays (vet_sc_set = gd: the vet_gd set, pole on x1), lagged per hydro step, FULL tensor D = K/J (eigenvalue floor 0.000000), vet_bc_diff = face`
  - `Terminating on cycle limit`
### G4_B1: s/cycle nan (cycles ); ; 
- elapsed: c0 0.469
  - `<rad_vet> vet_gd direction set: 160 directions; moment errors |sum w - 1| 4.440892e-16, |sum w n| 2.612927e-17, |sum w nn - I/3| 1.665335e-16 (rotated 3.885781e-16), 4th-moment anisotropy 2.775558e-16 (rotated 1.387779e-16)`
  - `### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel6/src/rad_vet/rays/sc_halo.hpp`
  - `### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel6/src/rad_vet/rays/sc_halo.hpp`
### G4_A2: s/cycle 0.2128 (cycles 2-6); cpu time used  = 1.659391e+00; zone-cycles/cpu_second = 7.582852e+06
- elapsed: c0 0.499; c1 1.091; c2 1.307; c3 1.520; c4 1.739; c5 1.951; c6 2.158
  - `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=9.774438e-01 (1.396348e-01 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=9.774438e+299`
  - `<rad_vet> vet_sc_set = gd: 160 rays, min |mu1| 0.183435, x1 mirror pairs yes`
  - `<rad_vet> closure = sc: short characteristics, 160 rays (vet_sc_set = gd: the vet_gd set, pole on x1), lagged per hydro step, FULL tensor D = K/J (eigenvalue floor 0.000000), vet_bc_diff = face`
  - `Terminating on cycle limit`
### G4_B2: s/cycle nan (cycles ); ; 
- elapsed: c0 0.471
  - `<rad_vet> vet_gd direction set: 160 directions; moment errors |sum w - 1| 4.440892e-16, |sum w n| 2.612927e-17, |sum w nn - I/3| 1.665335e-16 (rotated 3.885781e-16), 4th-moment anisotropy 2.775558e-16 (rotated 1.387779e-16)`
  - `### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel6/src/rad_vet/rays/sc_halo.hpp`
  - `### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel6/src/rad_vet/rays/sc_halo.hpp`
### G4x_A: s/cycle 0.1633 (cycles 2-6); cpu time used  = 1.453649e+00; zone-cycles/cpu_second = 8.656088e+06
- elapsed: c0 0.554; c1 1.188; c2 1.354; c3 1.518; c4 1.680; c5 1.850; c6 2.007
  - `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=6.383429e-01 (9.119184e-02 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=6.383429e+299`
  - `<rad_vet> vet_sc_set = gd: 160 rays, min |mu1| 0.183435, x1 mirror pairs yes`
  - `<rad_vet> closure = sc: short characteristics, 160 rays (vet_sc_set = gd: the vet_gd set, pole on x1), lagged per hydro step, FULL tensor D = K/J (eigenvalue floor 0.000000), vet_bc_diff = face`
  - `Terminating on cycle limit`
### G4x_B: s/cycle nan (cycles ); ; 
- elapsed: 
  - `### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel6/src/rad_vet/rad_vet_geometry.cpp at line 362`
  - `### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel6/src/rad_vet/rad_vet_geometry.cpp at line 362`
### G4L_A: s/cycle 0.6517 (cycles 2-6); cpu time used  = 4.557785e+00; zone-cycles/cpu_second = 2.208601e+07
- elapsed: c0 0.522; c1 1.790; c2 2.472; c3 3.138; c4 3.798; c5 4.451; c6 5.079
  - `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=2.052772e+00 (2.932531e-01 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=2.052772e+300`
  - `<rad_vet> vet_sc_set = gd: 160 rays, min |mu1| 0.183435, x1 mirror pairs yes`
  - `<rad_vet> closure = sc: short characteristics, 160 rays (vet_sc_set = gd: the vet_gd set, pole on x1), lagged per hydro step, FULL tensor D = K/J (eigenvalue floor 0.000000), vet_bc_diff = face`
  - `Terminating on cycle limit`
### G4L_B: s/cycle nan (cycles ); ; 
- elapsed: c0 0.495
  - `<rad_vet> vet_gd direction set: 160 directions; moment errors |sum w - 1| 4.440892e-16, |sum w n| 2.612927e-17, |sum w nn - I/3| 1.665335e-16 (rotated 3.885781e-16), 4th-moment anisotropy 2.775558e-16 (rotated 1.387779e-16)`
  - `### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel6/src/rad_vet/rays/sc_halo.hpp`
  - `### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel6/src/rad_vet/rays/sc_halo.hpp`
### kp_G1_A: s/cycle 0.3566 (cycles 2-6); cpu time used  = 2.382246e+00; zone-cycles/cpu_second = 5.281953e+06
- elapsed: c0 0.030; c1 0.623; c2 0.986; c3 1.349; c4 1.699; c5 2.062; c6 2.412
  - `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=1.095108e+00 (1.564440e-01 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=1.095108e+300`
  - `<rad_vet> vet_sc_set = gd: 160 rays, min |mu1| 0.183435, x1 mirror pairs yes`
  - `<rad_vet> closure = sc: short characteristics, 160 rays (vet_sc_set = gd: the vet_gd set, pole on x1), lagged per hydro step, FULL tensor D = K/J (eigenvalue floor 0.000000), vet_bc_diff = face`
  - `Terminating on cycle limit`
### kp_G1_B: s/cycle 0.5278 (cycles 2-6); cpu time used  = 7.136497e+00; zone-cycles/cpu_second = 1.763178e+06
- elapsed: c0 0.035; c1 4.512; c2 5.060; c3 5.594; c4 6.128; c5 6.658; c6 7.171
  - `<rad_vet> vet_gd: rank 0 fenced seconds: source 9.313888e-03, sweeps 1.923346e+00 (of which per-shell halo+walls 3.930890e-04, of which MPI wait 0.000000e+00), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
  - `<rad_vet> vet_gd halo (rank 0): exchanges 0.000000e+00, bytes sent 0.000000e+00 (0.000000e+00 per exchange), post 0.000000e+00 s, partners 0; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
  - `<rad_vet> vet_gd direction set: 160 directions; moment errors |sum w - 1| 4.440892e-16, |sum w n| 2.612927e-17, |sum w nn - I/3| 1.665335e-16 (rotated 3.885781e-16), 4th-moment anisotropy 2.775558e-16 (rotated 1.387779e-16)`
  - `<rad_vet> vet_gd: global-direction SC closure on the sp wedge, product GL8x20 (160 directions); z-rotation every 0 cycle(s) (0 = fixed), first angle 0.000000e+00; rank 0 wall ghost columns 62016, direction maps exact 3112960 of 9922560 (the rest: nearest direction)`
  - `Terminating on cycle limit`
### kp_G4_A: s/cycle 0.2375 (cycles 2-6); cpu time used  = 1.821074e+00; zone-cycles/cpu_second = 6.909610e+06
- elapsed: c0 0.500; c1 1.130; c2 1.371; c3 1.613; c4 1.853; c5 2.091; c6 2.321
  - `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=1.003778e+00 (1.433969e-01 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=1.003778e+300`
  - `<rad_vet> vet_sc_set = gd: 160 rays, min |mu1| 0.183435, x1 mirror pairs yes`
  - `<rad_vet> closure = sc: short characteristics, 160 rays (vet_sc_set = gd: the vet_gd set, pole on x1), lagged per hydro step, FULL tensor D = K/J (eigenvalue floor 0.000000), vet_bc_diff = face`
  - `Terminating on cycle limit`
### kp_G4_B: s/cycle nan (cycles ); ; 
- elapsed: c0 0.464
  - `<rad_vet> vet_gd direction set: 160 directions; moment errors |sum w - 1| 4.440892e-16, |sum w n| 2.612927e-17, |sum w nn - I/3| 1.665335e-16 (rotated 3.885781e-16), 4th-moment anisotropy 2.775558e-16 (rotated 1.387779e-16)`
  - `### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel6/src/rad_vet/rays/sc_halo.hpp`
  - `### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel6/src/rad_vet/rays/sc_halo.hpp`
- kp_G1_A kp.0 top 25:
```
# total 2.428666 s
      0.615664        896 vet_vet_sweep_ray
      0.464696        896 vet_vet_sweep_sum
      0.134942        141 cln_rhsaf
      0.125951       2130 cln_pcr2
      0.059981         94 cln_fsch
      0.044085        141 cln_fs0
      0.034663        280 cln_apply_st
      0.030740         48 cln_facecoef
      0.028513       1886 <deep_copy>
      0.028318        426 cls_restrict
      0.028220         47 cln_rhs
      0.026974         94 cln_fs2
      0.026716         36 cln_velocity
      0.023918        493 vet_lay_red2
      0.023239        123 cln_cslat
      0.021949        108 cln_fwc
      0.021862        141 cln_fwb
      0.021635         12 dfl_th
      0.020836       1362 cln_ghost1
      0.019933         47 cln_drag
      0.019286       1374 cln_hd_copy
      0.017385         99 cln_dot_x
      0.016156        282 cln_fwbs
      0.015686         33 cln_wbxc
```
- kp_G1_B kp.0 top 25:
```
# total 7.293190 s
      1.905542       1792 vet_vgd_shell
      0.155676       2620 cln_pcr2
      0.155222        144 cln_rhsaf
      0.069745          7 cln_gd_fs
      0.068976        144 cln_fs0
      0.040962        331 cln_apply_st
      0.040368       2202 <deep_copy>
      0.035819        524 cls_restrict
      0.034077         37 cln_velocity
      0.034026         48 cln_fsch
      0.033025         48 cln_rhs
      0.030465         48 cln_facecoef
      0.029266        124 cln_cslat
      0.028814        604 vet_lay_red2
      0.028603         48 cln_drag
      0.028110          7 cln_gd
      0.028088        111 cln_fwc
      0.027488         96 cln_fs2
      0.026522         12 dfl_th
      0.024784       1515 cln_ghost1
      0.023219        144 cln_fwb
      0.021572        123 cln_dot_x
      0.019717        111 cln_xdonor
      0.018976         37 cln_stencil
```
- kp_G4_A kp.0 top 25:
```
# total 2.360754 s
      0.886525        896 vet_vet_mb_sweep
      0.059260       2070 cln_pcr2
      0.039451        144 cln_rhsaf
      0.030683        192 cln_fsch
      0.022794        481 vet_lay_red2
      0.018999       1922 <deep_copy>
      0.016872       1482 cln_hd_pack
      0.016752       1482 cln_hd_unpack
      0.015807       1470 cln_ghost1
      0.015178        144 cln_fs0
      0.012744         96 cln_fs2
      0.011879        414 cls_restrict
      0.010967        274 cln_apply_st
      0.010759        889 vet_vet_plane_unpka
      0.010600        889 vet_vet_plane_packa
      0.010292        889 vet_vet_plane_loc
      0.010270         12 hflux_x1
      0.008783         48 cln_facecoef
      0.008570         33 cln_wbxc
      0.008281         48 cln_rhs
      0.007948        144 cln_fwb
      0.007702        111 cln_fwc
      0.007540         37 cln_velocity
      0.007119         11 cln_wbprb
```

## scspeed summary and an extra layout (02-17, same binary)
- The sp arms of G4, G4x and G4L FATAL at start-up by design, so the note's lateral split cannot run on sp:
  - G4 / G4L (2 x 2 lateral blocks): `vet_gd_theta_wall = mirror needs ONE MeshBlock along x2 (theta)`.
  - G4x (2 blocks along x1): `on a spherical-polar mesh supports only ... one MeshBlock along x1`.
  - The Cartesian arms of all three ran.
- Extra arms, not in the note (cheap, run so the 4-GPU comparison exists): **G4p** = 128^3 split into 4 blocks along phi only (meshblock 128 x 128 x 32), order A B A B; **G4Lp** = 256^3 with meshblock 256 x 256 x 64, A then B. Everything else as in the note.
- s/cycle over cycles 2-6 and the timer lines (rank 0):
| run | s/cycle | Cartesian vet_sc | sp vet_gd fenced | sp halo | FATAL/CW |
|---|---|---|---|---|---|
| G1_A1 | 0.3358 | SC 1.076 s (0.1538 per call) |  |  | 0/0 |
| G1_B1 | 0.5027 |  | sweeps 1.908 s, source 0.0089, halo+walls 1.81e-04, MPI wait 0.00e+00 | exchanges 0 | 0/0 |
| G1_A2 | 0.3362 | SC 1.076 s (0.1537 per call) |  |  | 0/0 |
| G1_B2 | 0.5028 |  | sweeps 1.907 s, source 0.0088, halo+walls 1.68e-04, MPI wait 0.00e+00 | exchanges 0 | 0/0 |
| G4_A1 | 0.2156 | SC 0.990 s (0.1414 per call) |  |  | 0/0 |
| G4_A2 | 0.2128 | SC 0.977 s (0.1396 per call) |  |  | 0/0 |
| G4x_A | 0.1633 | SC 0.638 s (0.0912 per call) |  |  | 0/0 |
| G4L_A | 0.6517 | SC 2.053 s (0.2933 per call) |  |  | 0/0 |
| G4p_A1 | 0.2184 | SC 0.989 s (0.1413 per call) |  |  | 0/0 |
| G4p_B1 | 0.2436 |  | sweeps 2.697 s, source 0.2119, halo+walls 2.67e+00, MPI wait 3.12e-01 | exchanges 1 | 0/0 |
| G4p_A2 | 0.2184 | SC 0.989 s (0.1412 per call) |  |  | 0/0 |
| G4p_B2 | 0.2432 |  | sweeps 2.681 s, source 0.1655, halo+walls 2.66e+00, MPI wait 3.13e-01 | exchanges 1 | 0/0 |
| G4Lp_A | 0.6693 | SC 2.083 s (0.2975 per call) |  |  | 0/0 |
| G4Lp_B | 1.1626 |  | sweeps 11.281 s, source 0.2040, halo+walls 1.12e+01, MPI wait 7.55e-01 | exchanges 3 | 0/0 |
- Reading:
  - 1 GPU, 128^3: the Cartesian sweep is faster per cycle, 0.336 against 0.503 s/cycle (sp 1.50x). The timers agree: SC 1.08 s against sp sweeps 1.91 s over the run (1.77x).
  - 4 GPUs, phi split: 0.218 against 0.243 s/cycle (sp 1.12x). The sp sweep time is almost all per-shell halo+walls (2.66 of 2.68 s); MPI wait is only 0.31 s.
  - 4 GPUs, 256^3 phi split: 0.669 against 1.163 s/cycle (sp 1.74x). sp sweeps take 11.3 s, of which halo+walls 11.2 s.
  - So on more than one GPU the curved path's per-shell halo dominates its cost.
- kp tool: `kp_G1_A/B` and `kp_G4_A` kp.0 were posted above (kp_G4_B FATAL, as G4_B). No rst/bin written.

## Full FATAL text of G4_B1 (your 27zi): `grep -A3 FATAL runs/G4_B1/log`; all 4 ranks print the same 2 lines, nothing else follows
```
### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel6/src/rad_vet/rays/sc_halo.hpp
<rad_vet>/vet_gd: vet_gd_theta_wall = mirror needs ONE MeshBlock along x2 (theta)
### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel6/src/rad_vet/rays/sc_halo.hpp
<rad_vet>/vet_gd: vet_gd_theta_wall = mirror needs ONE MeshBlock along x2 (theta)
```
- G4L_B and kp_G4_B print the same. The message is complete: it is the theta mirror wall (reflecting theta boundaries), which needs one MeshBlock along x2, as you guessed. Round 3 (sp2/cart2) follows after the LZAB arms.

## scspeed round 3 (your 27zi): sp2 / cart2 decks at 66b1c4ef4, node 02-11, suite binary rel6 8cc98fb3, one rank per GPU
```
S1 rc=0 wall=11 s args: sp2.athinput time/nlim=1
S4 rc=0 wall=9 s args: sp2.athinput time/nlim=1 meshblock/nx2=32 meshblock/nx3=32
H1_A1 rc=0 wall=7 s args: cart2.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=128 meshblock/nx3=128 rad_vet/kappa_p=32
H1_B1 rc=0 wall=7 s args: sp2.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=128 meshblock/nx3=128 rad_vet/kappa_p=11.851851851851851
H1_A2 rc=0 wall=12 s args: cart2.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=128 meshblock/nx3=128 rad_vet/kappa_p=32
H1_B2 rc=0 wall=13 s args: sp2.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=128 meshblock/nx3=128 rad_vet/kappa_p=11.851851851851851
H4_A1 rc=0 wall=10 s args: cart2.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=64 meshblock/nx3=64 rad_vet/kappa_p=32
H4_B1 rc=0 wall=16 s args: sp2.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=64 meshblock/nx3=64 rad_vet/kappa_p=11.851851851851851
H4_A2 rc=0 wall=16 s args: cart2.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=64 meshblock/nx3=64 rad_vet/kappa_p=32
H4_B2 rc=0 wall=10 s args: sp2.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=64 meshblock/nx3=64 rad_vet/kappa_p=11.851851851851851
H4L_A rc=0 wall=18 s args: cart2.athinput mesh/nx1=256 mesh/nx2=256 mesh/nx3=256 meshblock/nx1=256 meshblock/nx2=128 meshblock/nx3=128 rad_vet/kappa_p=64
H4L_B rc=0 wall=19 s args: sp2.athinput mesh/nx1=256 mesh/nx2=256 mesh/nx3=256 meshblock/nx1=256 meshblock/nx2=128 meshblock/nx3=128 rad_vet/kappa_p=23.703703703703702
kp_H1_B rc=0 wall=7 s args: sp2.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=128 meshblock/nx3=128 rad_vet/kappa_p=11.851851851851851
kp_H4_A rc=0 wall=16 s args: cart2.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=64 meshblock/nx3=64 rad_vet/kappa_p=32
kp_H4_B rc=0 wall=11 s args: sp2.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=64 meshblock/nx3=64 rad_vet/kappa_p=11.851851851851851
written rst/bin: 0
```
### S1: s/cycle nan over cycles 2-6; cpu time used  = 5.948387e-01; zone-cycles/cpu_second = 4.406976e+05
- elapsed: c0 0.197; c1 0.792
  - `<rad_vet> vet_gd: rank 0 fenced seconds: source 2.209747e-02, sweeps 3.946630e-02 (of which per-shell halo+walls 3.027000e-05, of which MPI wait 0.000000e+00), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
  - `<rad_vet> vet_gd halo (rank 0): exchanges 0.000000e+00, bytes sent 0.000000e+00 (0.000000e+00 per exchange), post 0.000000e+00 s, partners 0; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
  - `<rad_vet> vet_gd direction set: 126 directions; moment errors |sum w - 1| 6.661338e-16, |sum w n| 1.994932e-17, |sum w nn - I/3| 2.775558e-16 (rotated 1.110223e-16), 4th-moment anisotropy 2.498002e-16 (rotated 1.110223e-16)`
  - `<rad_vet> vet_gd: global-direction SC closure on the sp wedge, product GL7x18 staggered (126 directions); z-rotation every 0 cycle(s) (0 = fixed), first angle 0.000000e+00; rank 0 wall ghost columns 14400, direction maps exact 580608 of 1814400 (the rest: nearest direction)`
  - `Terminating on cycle limit`
### S4: s/cycle nan over cycles 2-6; cpu time used  = 7.091144e-01; zone-cycles/cpu_second = 3.696780e+05
- elapsed: c0 0.504; c1 1.213
  - `<rad_vet> vet_gd: rank 0 fenced seconds: source 5.116580e-02, sweeps 2.783647e-01 (of which per-shell halo+walls 2.570163e-01, of which MPI wait 7.622265e-03), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
  - `<rad_vet> vet_gd halo (rank 0): exchanges 2.600000e+02, bytes sent 4.797966e+08 (1.845372e+06 per exchange), post 0.000000e+00 s, partners 3; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
  - `<rad_vet> vet_gd direction set: 126 directions; moment errors |sum w - 1| 6.661338e-16, |sum w n| 1.994932e-17, |sum w nn - I/3| 2.775558e-16 (rotated 1.110223e-16), 4th-moment anisotropy 2.498002e-16 (rotated 1.110223e-16)`
  - `<rad_vet> vet_gd: global-direction SC closure on the sp wedge, product GL7x18 staggered (126 directions); z-rotation every 0 cycle(s) (0 = fixed), first angle 0.000000e+00; rank 0 wall ghost columns 6720, direction maps exact 322560 of 846720 (the rest: nearest direction)`
  - `Terminating on cycle limit`
### H1_A1: s/cycle 0.3063 over cycles 2-6; cpu time used  = 2.057489e+00; zone-cycles/cpu_second = 6.115664e+06
- elapsed: c0 0.031; c1 0.543; c2 0.863; c3 1.180; c4 1.482; c5 1.786; c6 2.088
  - `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=8.905608e-01 (1.272230e-01 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=8.905608e+299`
  - `<rad_vet> vet_sc_set = gd: 126 rays, min |mu1| 0.0851519, x1 mirror pairs no`
  - `Terminating on cycle limit`
### H1_B1: s/cycle 0.2867 over cycles 2-6; cpu time used  = 2.246501e+00; zone-cycles/cpu_second = 5.601116e+06
- elapsed: c0 0.039; c1 0.840; c2 1.139; c3 1.429; c4 1.715; c5 2.000; c6 2.286
  - `<rad_vet> vet_gd: rank 0 fenced seconds: source 6.601630e-03, sweeps 4.404852e-01 (of which per-shell halo+walls 1.660390e-04, of which MPI wait 0.000000e+00), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
  - `<rad_vet> vet_gd halo (rank 0): exchanges 0.000000e+00, bytes sent 0.000000e+00 (0.000000e+00 per exchange), post 0.000000e+00 s, partners 0; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
  - `<rad_vet> vet_gd direction set: 126 directions; moment errors |sum w - 1| 6.661338e-16, |sum w n| 1.994932e-17, |sum w nn - I/3| 2.775558e-16 (rotated 1.110223e-16), 4th-moment anisotropy 2.498002e-16 (rotated 1.110223e-16)`
  - `<rad_vet> vet_gd: global-direction SC closure on the sp wedge, product GL7x18 staggered (126 directions); z-rotation every 0 cycle(s) (0 = fixed), first angle 0.000000e+00; rank 0 wall ghost columns 35600, direction maps exact 1612800 of 4485600 (the rest: nearest direction)`
  - `Terminating on cycle limit`
### H1_A2: s/cycle 0.3065 over cycles 2-6; cpu time used  = 2.051783e+00; zone-cycles/cpu_second = 6.132671e+06
- elapsed: c0 0.030; c1 0.539; c2 0.856; c3 1.173; c4 1.475; c5 1.780; c6 2.082
  - `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=8.873419e-01 (1.267631e-01 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=8.873419e+299`
  - `<rad_vet> vet_sc_set = gd: 126 rays, min |mu1| 0.0851519, x1 mirror pairs no`
  - `Terminating on cycle limit`
### H1_B2: s/cycle 0.2856 over cycles 2-6; cpu time used  = 2.240059e+00; zone-cycles/cpu_second = 5.617224e+06
- elapsed: c0 0.033; c1 0.833; c2 1.131; c3 1.419; c4 1.705; c5 1.988; c6 2.273
  - `<rad_vet> vet_gd: rank 0 fenced seconds: source 6.710301e-03, sweeps 4.412851e-01 (of which per-shell halo+walls 1.767580e-04, of which MPI wait 0.000000e+00), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
  - `<rad_vet> vet_gd halo (rank 0): exchanges 0.000000e+00, bytes sent 0.000000e+00 (0.000000e+00 per exchange), post 0.000000e+00 s, partners 0; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
  - `<rad_vet> vet_gd direction set: 126 directions; moment errors |sum w - 1| 6.661338e-16, |sum w n| 1.994932e-17, |sum w nn - I/3| 2.775558e-16 (rotated 1.110223e-16), 4th-moment anisotropy 2.498002e-16 (rotated 1.110223e-16)`
  - `<rad_vet> vet_gd: global-direction SC closure on the sp wedge, product GL7x18 staggered (126 directions); z-rotation every 0 cycle(s) (0 = fixed), first angle 0.000000e+00; rank 0 wall ghost columns 35600, direction maps exact 1612800 of 4485600 (the rest: nearest direction)`
  - `Terminating on cycle limit`
### H4_A1: s/cycle 0.1882 over cycles 2-6; cpu time used  = 1.503626e+00; zone-cycles/cpu_second = 8.368382e+06
- elapsed: c0 0.496; c1 1.054; c2 1.247; c3 1.437; c4 1.625; c5 1.815; c6 1.999
  - `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=8.172215e-01 (1.167459e-01 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=8.172215e+299`
  - `<rad_vet> vet_sc_set = gd: 126 rays, min |mu1| 0.0851519, x1 mirror pairs no`
  - `Terminating on cycle limit`
### H4_B1: s/cycle 0.1161 over cycles 2-6; cpu time used  = 1.362632e+00; zone-cycles/cpu_second = 9.234267e+06
- elapsed: c0 0.495; c1 1.263; c2 1.393; c3 1.513; c4 1.629; c5 1.743; c6 1.858
  - `<rad_vet> vet_gd: rank 0 fenced seconds: source 6.716173e-02, sweeps 2.682702e-01 (of which per-shell halo+walls 1.962171e-02, of which MPI wait 1.141048e-02), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
  - `<rad_vet> vet_gd halo (rank 0): exchanges 1.813000e+03, bytes sent 6.990299e+08 (3.855653e+05 per exchange), post 7.327997e-03 s, partners 3; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
  - `<rad_vet> vet_gd direction set: 126 directions; moment errors |sum w - 1| 6.661338e-16, |sum w n| 1.994932e-17, |sum w nn - I/3| 2.775558e-16 (rotated 1.110223e-16), 4th-moment anisotropy 2.498002e-16 (rotated 1.110223e-16)`
  - `<rad_vet> vet_gd: global-direction SC closure on the sp wedge, product GL7x18 staggered (126 directions); z-rotation every 0 cycle(s) (0 = fixed), first angle 0.000000e+00; rank 0 wall ghost columns 13900, direction maps exact 718200 of 1751400 (the rest: nearest direction)`
  - `Terminating on cycle limit`
### H4_A2: s/cycle 0.1869 over cycles 2-6; cpu time used  = 1.485051e+00; zone-cycles/cpu_second = 8.473048e+06
- elapsed: c0 0.488; c1 1.034; c2 1.225; c3 1.416; c4 1.600; c5 1.790; c6 1.973
  - `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=8.119568e-01 (1.159938e-01 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=8.119568e+299`
  - `<rad_vet> vet_sc_set = gd: 126 rays, min |mu1| 0.0851519, x1 mirror pairs no`
  - `Terminating on cycle limit`
### H4_B2: s/cycle 0.1153 over cycles 2-6; cpu time used  = 1.372168e+00; zone-cycles/cpu_second = 9.170097e+06
- elapsed: c0 0.488; c1 1.271; c2 1.400; c3 1.519; c4 1.634; c5 1.748; c6 1.861
  - `<rad_vet> vet_gd: rank 0 fenced seconds: source 6.307289e-02, sweeps 2.795505e-01 (of which per-shell halo+walls 1.936070e-02, of which MPI wait 1.134798e-02), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
  - `<rad_vet> vet_gd halo (rank 0): exchanges 1.813000e+03, bytes sent 6.990299e+08 (3.855653e+05 per exchange), post 7.767704e-03 s, partners 3; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
  - `<rad_vet> vet_gd direction set: 126 directions; moment errors |sum w - 1| 6.661338e-16, |sum w n| 1.994932e-17, |sum w nn - I/3| 2.775558e-16 (rotated 1.110223e-16), 4th-moment anisotropy 2.498002e-16 (rotated 1.110223e-16)`
  - `<rad_vet> vet_gd: global-direction SC closure on the sp wedge, product GL7x18 staggered (126 directions); z-rotation every 0 cycle(s) (0 = fixed), first angle 0.000000e+00; rank 0 wall ghost columns 13900, direction maps exact 718200 of 1751400 (the rest: nearest direction)`
  - `Terminating on cycle limit`
### H4L_A: s/cycle 0.6065 over cycles 2-6; cpu time used  = 4.238243e+00; zone-cycles/cpu_second = 2.375119e+07
- elapsed: c0 0.520; c1 1.696; c2 2.332; c3 2.950; c4 3.563; c5 4.178; c6 4.758
  - `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=1.735953e+00 (2.479933e-01 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=1.735953e+300`
  - `<rad_vet> vet_sc_set = gd: 126 rays, min |mu1| 0.0851519, x1 mirror pairs no`
  - `Terminating on cycle limit`
### H4L_B: s/cycle 0.5943 over cycles 2-6; cpu time used  = 5.036487e+00; zone-cycles/cpu_second = 1.998681e+07
- elapsed: c0 0.510; c1 2.533; c2 3.169; c3 3.775; c4 4.371; c5 4.958; c6 5.547
  - `<rad_vet> vet_gd: rank 0 fenced seconds: source 8.540212e-02, sweeps 1.268784e+00 (of which per-shell halo+walls 1.777907e-01, of which MPI wait 2.456280e-02), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
  - `<rad_vet> vet_gd halo (rank 0): exchanges 3.619000e+03, bytes sent 3.168571e+09 (8.755377e+05 per exchange), post 1.572514e-02 s, partners 3; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
  - `<rad_vet> vet_gd direction set: 126 directions; moment errors |sum w - 1| 6.661338e-16, |sum w n| 1.994932e-17, |sum w nn - I/3| 2.775558e-16 (rotated 1.110223e-16), 4th-moment anisotropy 2.498002e-16 (rotated 1.110223e-16)`
  - `<rad_vet> vet_gd: global-direction SC closure on the sp wedge, product GL7x18 staggered (126 directions); z-rotation every 0 cycle(s) (0 = fixed), first angle 0.000000e+00; rank 0 wall ghost columns 31947, direction maps exact 1712718 of 4025322 (the rest: nearest direction)`
  - `Terminating on cycle limit`
### kp_H1_B: s/cycle 0.3072 over cycles 2-6; cpu time used  = 2.392797e+00; zone-cycles/cpu_second = 5.258662e+06
- elapsed: c0 0.035; c1 0.879; c2 1.199; c3 1.510; c4 1.817; c5 2.120; c6 2.428
  - `<rad_vet> vet_gd: rank 0 fenced seconds: source 6.820108e-03, sweeps 4.560092e-01 (of which per-shell halo+walls 3.330880e-04, of which MPI wait 0.000000e+00), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
  - `<rad_vet> vet_gd halo (rank 0): exchanges 0.000000e+00, bytes sent 0.000000e+00 (0.000000e+00 per exchange), post 0.000000e+00 s, partners 0; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
  - `<rad_vet> vet_gd direction set: 126 directions; moment errors |sum w - 1| 6.661338e-16, |sum w n| 1.994932e-17, |sum w nn - I/3| 2.775558e-16 (rotated 1.110223e-16), 4th-moment anisotropy 2.498002e-16 (rotated 1.110223e-16)`
  - `<rad_vet> vet_gd: global-direction SC closure on the sp wedge, product GL7x18 staggered (126 directions); z-rotation every 0 cycle(s) (0 = fixed), first angle 0.000000e+00; rank 0 wall ghost columns 35600, direction maps exact 1612800 of 4485600 (the rest: nearest direction)`
  - `Terminating on cycle limit`
### kp_H4_A: s/cycle 0.2120 over cycles 2-6; cpu time used  = 1.657317e+00; zone-cycles/cpu_second = 7.592338e+06
- elapsed: c0 0.492; c1 1.083; c2 1.301; c3 1.518; c4 1.727; c5 1.942; c6 2.149
  - `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=8.430594e-01 (1.204371e-01 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=8.430594e+299`
  - `<rad_vet> vet_sc_set = gd: 126 rays, min |mu1| 0.0851519, x1 mirror pairs no`
  - `Terminating on cycle limit`
### kp_H4_B: s/cycle 0.1470 over cycles 2-6; cpu time used  = 1.598915e+00; zone-cycles/cpu_second = 7.869654e+06
- elapsed: c0 0.501; c1 1.355; c2 1.513; c3 1.665; c4 1.812; c5 1.956; c6 2.100
  - `<rad_vet> vet_gd: rank 0 fenced seconds: source 6.568841e-02, sweeps 3.624421e-01 (of which per-shell halo+walls 1.123011e-01, of which MPI wait 1.211569e-02), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
  - `<rad_vet> vet_gd halo (rank 0): exchanges 1.813000e+03, bytes sent 6.990299e+08 (3.855653e+05 per exchange), post 8.705500e-03 s, partners 3; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
  - `<rad_vet> vet_gd direction set: 126 directions; moment errors |sum w - 1| 6.661338e-16, |sum w n| 1.994932e-17, |sum w nn - I/3| 2.775558e-16 (rotated 1.110223e-16), 4th-moment anisotropy 2.498002e-16 (rotated 1.110223e-16)`
  - `<rad_vet> vet_gd: global-direction SC closure on the sp wedge, product GL7x18 staggered (126 directions); z-rotation every 0 cycle(s) (0 = fixed), first angle 0.000000e+00; rank 0 wall ghost columns 13900, direction maps exact 718200 of 1751400 (the rest: nearest direction)`
  - `Terminating on cycle limit`
- kp_H1_B kp.0 top 25:
```
# total 2.503617 s
      0.439622       1792 vet_vgd_shell
      0.161634        147 cln_rhsaf
      0.112877       1910 cln_pcr2
      0.069569        147 cln_fs0
      0.068445         98 cln_fsch
      0.055398          7 cln_gd_fs
      0.053078       2110 <deep_copy>
      0.035015         38 cln_velocity
      0.033752         49 cln_rhs
      0.032356        261 cln_apply_st
      0.030450         48 cln_facecoef
      0.029463        125 cln_cslat
      0.029183         49 cln_drag
      0.029019        114 cln_fwc
      0.028016         98 cln_fs2
      0.026616         12 dfl_th
      0.026062        382 cls_restrict
      0.023926          7 cln_gd
      0.023826        147 cln_fwb
      0.023636       1459 cln_ghost1
      0.022636        463 vet_lay_red2
      0.021155       1471 cln_hd_copy
      0.020241        114 cln_xdonor
      0.019553         38 cln_stencil
```
- kp_H4_A kp.0 top 25:
```
# total 2.191280 s
      0.721682        896 vet_vet_mb_sweep
      0.058769       2070 cln_pcr2
      0.039397        141 cln_rhsaf
      0.031913        188 cln_fsch
      0.022810        480 vet_lay_red2
      0.018909       1901 <deep_copy>
      0.016641       1460 cln_hd_pack
      0.016419       1460 cln_hd_unpack
      0.015600       1448 cln_ghost1
      0.014962        141 cln_fs0
      0.012288         94 cln_fs2
      0.012007        889 vet_vet_plane_unpka
      0.011824        889 vet_vet_plane_packa
      0.011806        414 cls_restrict
      0.011275        889 vet_vet_plane_loc
      0.010914        273 cln_apply_st
      0.010511         12 hflux_x1
      0.008705         48 cln_facecoef
      0.008153         47 cln_rhs
      0.007793        141 cln_fwb
      0.007718         33 cln_wbxc
      0.007465        108 cln_fwc
      0.007344         36 cln_velocity
      0.006628        123 cln_cslat
```
- kp_H4_B kp.0 top 25:
```
# total 2.206511 s
      0.107492       1792 vet_vgd_shell
      0.061639       2090 cln_pcr2
      0.045703        147 cln_rhsaf
      0.036194       1820 vet_vgd_ring
      0.034643        196 cln_fsch
      0.033198       2486 <deep_copy>
      0.023973        499 vet_lay_red2
      0.021804        147 cln_fs0
      0.019674       1785 vet_vgd_hl_scatter
      0.019373       1785 vet_vgd_ipc_pack
      0.019181       1584 cln_hd_pack
      0.018947       1584 cln_hd_unpack
      0.018922       1792 vet_vgd_ipc_k2
      0.018002       1572 cln_ghost1
      0.013095         98 cln_fs2
      0.012842          7 cln_gd_fs
      0.012834        382 vet_vgd_hxb_mark
      0.012581        418 cls_restrict
      0.011311        279 cln_apply_st
      0.010569         12 hflux_x1
      0.010472          2 vet_vgd_shl
      0.009612        114 cln_fwc
      0.009578         49 cln_rhs
      0.009426         38 cln_velocity
```
