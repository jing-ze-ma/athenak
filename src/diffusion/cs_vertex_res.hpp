#ifndef DIFFUSION_CS_VERTEX_RES_HPP_
#define DIFFUSION_CS_VERTEX_RES_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cs_vertex_res.hpp
//! \brief cornerfix-1103 <mesh>/cs_vertex_3cell: the cube-vertex radial edge of the
//! cubed-sphere resistive operators (Ohmic, ambipolar; also for the seam current of
//! cs_seam_res = mean_j).  At a cube vertex only three cells meet: this block's corner
//! cell (kv, jv) and the two seam ghost cells (kv, jn) (across the x2 seam) and (kn, jv)
//! (across the x3 seam); the fourth cell of the Cartesian 4-cell stencil (kn, jn) is the
//! made-up corner block.  CsVertexAt says whether the radial edge (k, j) of block m is a
//! cube vertex and gives the three cells; CsVertex3 is a cell field at the vertex point.

#include "athena.hpp"
#include "mesh/mesh.hpp"

struct CsVtx3 {
  int jv, kv, jn, kn;   // corner cell (kv, jv); seam ghosts (kv, jn) and (kn, jv)
};

KOKKOS_INLINE_FUNCTION
bool CsVertexAt(const RegionSize &z, const int k, const int j, const int js, const int je,
                const int ks, const int ke, CsVtx3 &v) {
  const bool jl = (j == js && z.x2min <= -1.0 + 1.0e-12);
  const bool ju = (j == je+1 && z.x2max >= 1.0 - 1.0e-12);
  const bool kl = (k == ks && z.x3min <= -1.0 + 1.0e-12);
  const bool ku = (k == ke+1 && z.x3max >= 1.0 - 1.0e-12);
  if (!((jl || ju) && (kl || ku))) return false;
  v.jv = jl ? js : je;  v.jn = jl ? js-1 : je+1;
  v.kv = kl ? ks : ke;  v.kn = kl ? ks-1 : ke+1;
  return true;
}

//! cornerfix2-1103: a cell field AT the vertex point.  In this block's own (smoothly
//! extended) gnomonic chart the two seam ghost cells (kv, jn) and (kn, jv) are point
//! reflections of each other through the vertex, so their mean is the linear-exact,
//! convex (never overshooting) vertex value: second order, bounded at jumps.  (The
//! 3-cell mean with the corner cell is NOT linear-exact in this chart: first order.)
template <class A>
KOKKOS_INLINE_FUNCTION
Real CsVertex3(const A &a, const int m, const CsVtx3 &v, const int i) {
  return 0.5*(a(m,v.kv,v.jn,i) + a(m,v.kn,v.jv,i));
}
//! the same for a component of a 5-D array (bcc)
template <class A>
KOKKOS_INLINE_FUNCTION
Real CsVertex3n(const A &a, const int m, const int n, const CsVtx3 &v, const int i) {
  return 0.5*(a(m,n,v.kv,v.jn,i) + a(m,n,v.kn,v.jv,i));
}

#endif  // DIFFUSION_CS_VERTEX_RES_HPP_
