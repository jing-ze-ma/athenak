# NOTE 2026-11-27zzn orion -> Caltech: rayunify round 5, slab face scalars folded into the shell kernel (G1/G4 only)

From RAYUNIFY (orion). Thanks for 27zzm. With FS the cart2 rays cost 0.279 shell + 0.275 slab_phi + 0.178 slab_fs.

**FD** = branch **`caltech-rayunify-c9d5d3ee8`** = FS + one commit:
- the x2/x3 face scalars are computed in the slab shell kernel itself (no vet_vgd_slab_phi);
- the per-cell sums run once per pass (2 launches per build instead of 256).
- On Orion CPU all 605 payload files are bitwise vs FS, and nvcc compiles.
- Memory: the pass scratch is 5 x n_rays x cells x 8 B, i.e. 10.6 GB at 128^3 x 126 rays on one GPU. Please report peak GPU memory (nvidia-smi or the log's memory line) for G1 cart2.

Builds, built_in_pgens only:
- **FD**: default.
- **FD0**: the same commit with `-D CMAKE_CXX_FLAGS="-ffp-contract=off -DVET_SLAB_PHI_FOLD=0"` (separate build dir). Same bits; the phi kernel runs separately but the sums are still per pass. It isolates the cost of the fold.

## Arms
G1 (and G4 if time allows) of 27zzl: cart2 and sp2, 128^3, kp timer. Binaries FS, FD, FD0.

Report per run:
- the elapsed line;
- kp.0 top 10 (vet_vgd_shell, vet_vgd_slab_phi, vet_vgd_slab_fs with launch counts);
- for G1 cart2, the peak GPU memory.

-- RAYUNIFY (orion)
