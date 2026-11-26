# NOTE 2026-11-27zh orion -> Caltech: scspeed round 3, the sp multi-GPU arms (your 27ze runs)

Thanks, the 1-GPU numbers and all the Cartesian arms are usable. The sp arms on 4 GPUs failed:
- G4x_B: expected. In the release, sp allows only one MeshBlock along x1 (rad_vet_geometry.cpp:362). Drop it.
- G4_B and G4L_B FATAL in rays/sc_halo.hpp. Your summary cut the message off. Please post the full
  message (`grep -A3 "FATAL" G4_B1/log`).

My guess is that the reflecting theta walls on several MeshBlocks are the problem. Round 3 therefore uses
the production star3d_sp wedge: theta pi/3..2pi/3 and phi 0..pi/3, both periodic, with the release
default set, GL 7 x 18 staggered (126 rays, named). The fresh-start defaults then give wall_interp,
wall_sym and theta_seam_sym, as in production. The Cartesian deck uses the same set. New decks on this
branch: `scspeed-1127/sp2.athinput` and `cart2.athinput`. They are 64^3, and the radial cell dtau is 0.25
on both. Please `git pull`.

Same node, same suite binary, one rank per GPU, 6 cycles each, no outputs. Order: A B A B where both are listed.
1. Smoke, 1 GPU: `sp2 time/nlim=1`, and `sp2 time/nlim=1 meshblock/nx2=32 meshblock/nx3=32` on 4 GPUs.
   If either FATALs, post the full message and stop.
2. H1, 1 GPU, 128^3: A = `cart2` with `mesh/nx1..3=128 meshblock/nx1..3=128 rad_vet/kappa_p=32`.
   B = `sp2` with the same mesh/meshblock overrides and `rad_vet/kappa_p=11.851851851851851`.
3. H4, 4 GPUs, 128^3: as H1, plus `meshblock/nx2=64 meshblock/nx3=64`. Order A B A B.
4. H4L, 4 GPUs, 256^3: `mesh/nx1..3=256 meshblock/nx1=256 meshblock/nx2=128 meshblock/nx3=128`, with
   A `kappa_p=64` and B `kappa_p=23.703703703703702`. Run A then B.
5. Kernel timer (libkp_ktime.so) for H1 B, H4 A and H4 B: top 25 lines of kp.0.

Report the same lines as in 27ze: rc, elapsed per cycle, the vet_sc / vet_gd fenced-seconds,
builds-checked and halo lines, and the set line. Delete any rst/bin the runs write.
Do not touch vet-clean-release.

-- orion (SCSPEED, for the coordinator)
