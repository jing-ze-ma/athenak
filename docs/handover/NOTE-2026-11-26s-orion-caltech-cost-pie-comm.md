# NOTE 2026-11-26s orion -> Caltech: calculation against communication for the cost figure (addendum to 26r, the USER, 12:08 PST)

The figure the USER means is Fig. 20 of arXiv 2503.16627: an exploded donut, one coloured slice per module with a darker rim and its percentage, "cells (cores)" in the centre, and for the radiation module outer arcs that split it into calculations and communication (there: 19 / 17 / 11 %). The USER wants each module here (hydro, rays, rad-vet, others) also subdivided into **calculation and communication**. So in the profile of 26r please also return, per star and for the mean of the 4 ranks:
1. per module, GPU kernel time in two classes: calculation kernels, and pack / copy / unpack kernels (halo, boundary: CleanHaloPost / Finish "void", Pack / Unpack, corner fills, ray mark / shell exchange).
2. host-side MPI time per call site (Waitall, Isend, Test, Allreduce) by the site name the code prints (cln_hd_copy, cln_hd_unpack, vet_vgdcs_pack, vet_vgdcs_shell, cln_cslat, cpackx, SendBuff, ...): your mpit timer / gaps.txt; assign each site to hydro / rays / rad-vet / others, the unassigned one listed apart; GPU idle that is not MPI (sync, launch, memcpy) as its own "other" number.
3. the sum check: calculation + communication + idle = wall of the window (state the residual).

Use the per-cycle mean over the window, including the cycles with and without a ray rebuild. If GPU-resident copies cannot be split into calculation / communication with certainty, say which kernel names are ambiguous. Nothing else changes in 26r. The palette is still the USER's choice (no figures before that).

Also: please send me (and the USER) the URL of the palette page you staged.

-- orion (coordinator)
