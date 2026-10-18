#!/usr/bin/env python3
""" @package ./examples/idort/G25_setup/convert_ic.py
Convert the g25 cons12 IC (old IDORT-dev production, Fritz job 4325157) for the IDORT module (Config.sh here).

    python3 convert_ic.py IN.hdf5 OUT.hdf5 [--reduce] [--seed 42]

IN: rhoin-RSGg_golden25_p438_hr_bgr2_coreB_opal_c1s.ic.hdf5 (6.09M gas cells + one PartType1 point mass, 2.11e34 g,
at the box centre; 'Masses' holds the DENSITY, MESHRELAX_DENSITY_IN_INPUT; gas-only u; uniform NuclearComposition
H, He4, C12, N14, O16, Ne20 with X = 0.6538476353670762).

What is changed (everything else is copied as is):
  * BoundaryRole (BOUNDARY_CELLS), as I4_star_lite: the IC's background beyond ~9e14 cm is a regular 32^3 lattice
    (spacing BoxSize/32 = 5.47e14 cm, offset 1/2); its outermost lattice layer (index 0 or 31 in any dimension, 5768
    cells) becomes WALL_SOLID (frozen; vacuum for the radiation: IdortBndI is absent = 0), the next layer (index 1 or
    30, 5048 cells) WALL_FLUID (slip wall). These cells keep their background state (rho 1e-20 g/cm^3, u 4.84e10
    erg/g, T ~ 5e2 K: inside the EOS table, constant continuation of the opacity tables as for every CSM cell).
  * Header BoxSize (0 in the IC) is set to the cons12 BoxSize; the unused 'Density' dataset (overwritten by
    MESHRELAX_DENSITY_IN_INPUT) is dropped.
  * --reduce: smoke-test variant with ~1.1M cells. Cells are thinned at random per radial zone (keep fractions in
    REDUCE below; the 32^3 background lattice is always kept). Because the IC carries point values of rho and u of the
    MESA-derived profile at the generators, and the masses are rho V of the new mesh (MESHRELAX_DENSITY_IN_INPUT), a
    random subset is still a sampling of the same profile; only the resolution changes (zone spacing x keep^-1/3)
    and the mesh is less regular (Poisson thinning; the mesh regularisation relaxes it). The thinning keeps the
    shell structure of the photosphere layer (cells there are kept with the same fraction at every radius).
"""
import argparse

import h5py
import numpy as np

BOXSIZE = 1.75154570218571e+16  # cons12 param BoxSize [cm]
NLAT = 32                       # background lattice of the IC
ROLE_WALL_SOLID, ROLE_WALL_FLUID = 2, 4
# (r_max [cm], keep fraction): core and envelope / photosphere layer / CSM / outer background
REDUCE = [(5.70e13, 0.5), (5.87e13, 0.1), (3.0e14, 0.1), (5.0e14, 0.25), (np.inf, 1.0)]

ap = argparse.ArgumentParser()
ap.add_argument('inp')
ap.add_argument('out')
ap.add_argument('--reduce', action='store_true')
ap.add_argument('--seed', type=int, default=42)
a = ap.parse_args()

with h5py.File(a.inp, 'r') as f:
    g = f['PartType0']
    d0 = {k: g[k][:] for k in g if k != 'Density'}
    d1 = {k: f['PartType1'][k][:] for k in f['PartType1']}
    hattrs = dict(f['Header'].attrs)

pos = d0['Coordinates']
dx = BOXSIZE / NLAT
idx = pos / dx - 0.5
onlat = np.all(np.abs(idx - np.round(idx)) < 1e-6, axis=1)
ii = np.round(idx).astype(int)
edge = np.min(np.minimum(ii, NLAT - 1 - ii), axis=1)
role = np.zeros(len(pos), dtype=np.int32)
role[onlat & (edge == 0)] = ROLE_WALL_SOLID
role[onlat & (edge == 1)] = ROLE_WALL_FLUID
assert (role == ROLE_WALL_SOLID).sum() == NLAT**3 - (NLAT - 2)**3, 'outer lattice layer incomplete'
assert (role == ROLE_WALL_FLUID).sum() == (NLAT - 2)**3 - (NLAT - 4)**3, 'second lattice layer incomplete'
d0['BoundaryRole'] = role

if a.reduce:
    r = np.sqrt(np.sum((pos - 0.5 * BOXSIZE)**2, axis=1))
    rng = np.random.RandomState(a.seed)
    keep = np.zeros(len(pos), dtype=bool)
    rlo = 0.0
    for rmax, fk in REDUCE:
        z = (r >= rlo) & (r < rmax)
        keep |= z & (rng.uniform(size=len(pos)) < fk)
        print('zone %.3e - %.3e: %d -> %d cells' % (rlo, rmax, z.sum(), (z & keep).sum()))
        rlo = rmax
    keep |= onlat
    d0 = {k: v[keep] for k, v in d0.items()}

n0 = len(d0['ParticleIDs'])
with h5py.File(a.out, 'w') as f:
    h = f.create_group('Header')
    for k, v in hattrs.items():
        h.attrs[k] = v
    npart = np.array([n0, len(d1['ParticleIDs']), 0, 0, 0, 0], dtype=np.int32)
    h.attrs['NumPart_ThisFile'] = npart
    h.attrs['NumPart_Total'] = npart
    h.attrs['BoxSize'] = BOXSIZE
    for name, d in (('PartType0', d0), ('PartType1', d1)):
        g = f.create_group(name)
        for k, v in d.items():
            g[k] = v
print('%s: %d gas cells (%d WALL_SOLID, %d WALL_FLUID), %d point mass(es) %s g at %s' %
      (a.out, n0, (d0['BoundaryRole'] == 2).sum(), (d0['BoundaryRole'] == 4).sum(), len(d1['ParticleIDs']),
       d1['Masses'], d1['Coordinates']))
