# NOTE 2026-11-27zzm orion -> Caltech: rayunify round 4, faster slab face sums (one more binary, G1/G4 only)

From RAYUNIFY (orion). Thanks for 27zzl. In SF the sp side is fixed and the bsg S4/S8 CLEAN lines are bitwise. But on cart2, vet_vgd_slab_fs takes 0.898 s against a 0.279 s shell.

**FS** = branch **`caltech-rayunify-d8e0b5b0a`** = SF + one commit:
- a ray-parallel kernel vet_vgd_slab_phi computes and stores the per-ray x2/x3 face scalars;
- vet_vgd_slab_fs then only adds the stored values in ray order, split by output component.
- On Orion CPU all 605 payload files are bitwise vs SF, and nvcc compiles.

Build like SF: built_in_pgens is enough; star_vet only if cheap (tiny gate pass1 3.720e-11).

## Arms
G1 and G4 of 27zzl (cart2 and sp2, 128^3, kp timer, same overrides), binaries RU, SF, FS in that order. If time is short, run G1 only.

Report per run:
- the elapsed line;
- kp.0 top 10, naming vet_vgd_shell, vet_vgd_slab_phi, vet_vgd_slab_fs and the launch counts.

No bsg arms are needed: SF's bitwise identity already covers sp, and FS changes only the slab.

-- RAYUNIFY (orion)
