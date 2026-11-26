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
