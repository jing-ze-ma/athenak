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
#include "bvals/cs_seam_fc3.hpp"

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
//! cornerfix3-1103: raised to a THIRD-order point value by the pair's curvature term,
//! 1/8 of the second difference along the pair's direction centred on this block's own
//! corner cell (the cells (kn, jv+dj), (kv, jv), (kv+dk, jn): real, corner-free), gated
//! by the shared smoothness weight and bounded to the range of the five values
//! (cs_fc3::CurvCorrect): the pair mean at a jump, never an overshoot.
template <class F>
KOKKOS_INLINE_FUNCTION
Real CsVertexCell(const F &f, const CsVtx3 &v) {
  const int dj = (v.jn < v.jv) ? 1 : -1, dk = (v.kn < v.kv) ? 1 : -1;
  const Real p1 = f(v.kv, v.jn), p2 = f(v.kn, v.jv);
  const Real fa = f(v.kn, v.jv+dj), fb = f(v.kv, v.jv), fc = f(v.kv+dk, v.jn);
  const Real lo = fmin(fmin(fmin(p1, p2), fmin(fa, fb)), fc);
  const Real hi = fmax(fmax(fmax(p1, p2), fmax(fa, fb)), fc);
  return cs_fc3::CurvCorrect(0.5*(p1 + p2), fa, fb, fc, 0.125, lo, hi);
}
template <class A>
KOKKOS_INLINE_FUNCTION
Real CsVertex3(const A &a, const int m, const CsVtx3 &v, const int i) {
  return CsVertexCell([&](const int k, const int j) {return a(m,k,j,i);}, v);
}
//! the same for a component of a 5-D array (bcc)
template <class A>
KOKKOS_INLINE_FUNCTION
Real CsVertex3n(const A &a, const int m, const int n, const CsVtx3 &v, const int i) {
  return CsVertexCell([&](const int k, const int j) {return a(m,n,k,j,i);}, v);
}


//! cornerfix4-1103: a stencil read of cell (k, j) that may fall in the made-up cube-
//! vertex corner block (k outside [ks,ke] AND j outside [js,je] at a corner whose two
//! flanking block faces are panel edges): the parallelogram value
//!   f(k, jv) + f(kv, j) - f(kv, jv)   (jv, kv = the nearest active column / row)
//! from the seam ghost row cell, the seam ghost column cell and this block's own corner
//! cell -- linear-exact in the block's extended chart, so a difference stencil through
//! it stays second order (the same construction as the corner-cell edges of
//! CsVertexJFix).  Every other (k, j): f(k, j) unchanged.
KOKKOS_INLINE_FUNCTION
bool CsCornerRead(const RegionSize &z, const int k, const int j, const int js,
                  const int je, const int ks, const int ke, int &kv, int &jv) {
  const bool jl = (j < js), jh = (j > je), kl = (k < ks), kh = (k > ke);
  if (!((jl || jh) && (kl || kh))) return false;
  const bool pj = jl ? (z.x2min <= -1.0 + 1.0e-12) : (z.x2max >= 1.0 - 1.0e-12);
  const bool pk = kl ? (z.x3min <= -1.0 + 1.0e-12) : (z.x3max >= 1.0 - 1.0e-12);
  if (!(pj && pk)) return false;
  jv = jl ? js : je; kv = kl ? ks : ke;
  return true;
}
template <class F>
KOKKOS_INLINE_FUNCTION
Real CsCornerCell(const F &f, const RegionSize &z, const int k, const int j,
                  const int js, const int je, const int ks, const int ke) {
  int kv, jv;
  if (!CsCornerRead(z, k, j, js, je, ks, ke, kv, jv)) return f(k, j);
  return f(k, jv) + f(kv, j) - f(kv, jv);
}

#endif  // DIFFUSION_CS_VERTEX_RES_HPP_
