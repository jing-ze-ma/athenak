#!/bin/bash            # this line only there to enable syntax highlighting in this file

## examples/idort/G25_setup/Config.sh
## g25 RSG validation run with the IDORT module (comparison with the old IDORT-dev production cons12, Fritz job
## 4325157). Base: C6_star_piece_levels / I4_star_lite (VORONOI_CELL, moving mesh, BOUNDARY_CELLS vacuum layer at the
## box faces, SPHERICAL_SOURCE core heating, IDORT 24 angles + DSA), gravity as cons12 (self-gravity of the gas, one
## fixed point mass of type 1 with exact N^2 gravity). Classic-mesh variant: Config_classic.sh. See README.md.

#--------------------------------------- Basic operation mode of code
NTYPES=6
NSOFTTYPES=4                             # as cons12 (types 0, 1, 2 used; 3 for 3-5)

#--------------------------------------- Mesh: moving
VORONOI
VORONOI_CELL
SHAREDNODEMEMORY                         # required by VORONOI_CELL
VORONOI_DYNAMIC_UPDATE
REGULARIZE_MESH_CM_DRIFT
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED
REGULARIZE_MESH_FACE_ANGLE
BOUNDARY_CELLS                           # box faces: WALL_SOLID (radiation: vacuum) + WALL_FLUID layer (slip wall), as I4
MESHRELAX_DENSITY_IN_INPUT               # cons12 IC: 'Masses' holds the density (mass = rho V after the first mesh)

#--------------------------------------- Refinement (as cons12; star bands IdortRefine*, idort-g25ref)
# Remove this block for a fixed cell set (the 1.19M-cell smoke test ran without it, see README.md).
REFINEMENT_SPLIT_CELLS
REFINEMENT_MERGE_CELLS
REFINEMENT_VOLUME_LIMIT
REFINEMENT_AROUND_DM
REFINEMENT_AROUND_DM_NO_FORCED_DEREF     # as the old IDORT-dev (forced derefinement inside the softening disabled)

#--------------------------------------- Time integration: local time steps (as cons12)
# VORONOI_CELL with local steps needs TIMESTEPS_WAKEUP -> MHD + RIEMANN_HLLD (B = 0; cons12 also ran MHD + HLLD).
# Global steps: add FORCE_EQUAL_TIMESTEPS (and TIMESTEPS_WAKEUP, MHD, RIEMANN_HLLD may stay).
TIMESTEPS_WAKEUP
MHD
RIEMANN_HLLD
TREE_BASED_TIMESTEPS                     # as cons12 (radiation signal speed via the IDORT hook)

#--------------------------------------- Gravity: as cons12 (gas self-gravity + one fixed point mass, type 1)
SELFGRAVITY
HIERARCHICAL_GRAVITY
CELL_CENTER_GRAVITY
GRAVITY_NOT_PERIODIC
ALLOW_DIRECT_SUMMATION
DIRECT_SUMMATION_THRESHOLD=1000
EXACT_GRAVITY_FOR_PARTICLE_TYPE=1        # point mass: N^2 gravity
NO_SELFGRAVITY_TYPE=1
NO_GRAVITY_TYPE=1                        # the point mass feels no force ...
EXACT_GRAVITY_REACTION                   # ... but acts on the gas
GRAVITY_PM_ENERGY_CONSERVING             # conservative point-mass gravity energy (as run long_m8)
DMFIXED                                  # ... and is held at rest (velocity zeroed)
DMLOWESTTIMEBIN
EVALPOTENTIAL
MULTIPLE_NODE_SOFTENING
ADAPTIVE_HYDRO_SOFTENING
CHUNKING

#--------------------------------------- Matter: gas EOS from the module table (IdortEosMode 2), X per cell
EOS_PASSIVE                              # NuclearComposition carried; IDORT table EOS uses Composition[0] = X
EOS_NSPECIES=6                           # H, He4, C12, N14, O16, Ne20 (cons12 IC)
VARIABLE_GAMMA                           # IDORT EOS hook sets GammaE, GammaC

#--------------------------------------- Core heating (old IRT_STAR_CORE_OPTB heating source)
SPHERICAL_SOURCE

#--------------------------------------- IDORT
IDORT
IDORT_NUM_ANGLES=24                      # BRULS order 2 (cons12: IdortBrulsNum 2)
IDORT_DSA                                # recommended for star runs

#--------------------------------------- Precision and output
DOUBLEPRECISION=1
OUTPUT_IN_DOUBLEPRECISION
OUTPUT_VOLUME
OUTPUTPOTENTIAL
HAVE_HDF5
DEBUG
