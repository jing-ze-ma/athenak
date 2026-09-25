#ifndef RAD_M1_RAD_M1_PARFOR_HPP_
#define RAD_M1_RAD_M1_PARFOR_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_parfor.hpp
//! \brief par_for_lb: the 4-D athena par_for (same flattening, same index arithmetic)
//! with Kokkos::LaunchBounds<256,1> (m1-fast4, tests_m1/runs_5r_fast4).
//!
//! Kokkos' HIP backend compiles a plain RangePolicy kernel for up to 1024 threads per
//! block, which caps a kernel at 128 VGPRs on gfx942; the heavy per-cell kernels of the
//! implicit solve (EOS and opacity-table lookups inlined, 150-250 live registers) then
//! SPILL to scratch (m1_impl_i1 111 spilled VGPRs, m1_opacity 121, m1_vimp_j 69).  Kokkos
//! launches a RangePolicy with 256 threads per block anyway, so the bound only lets the
//! register allocator use what it needs.  The arithmetic is unchanged.

#include <string>

#include "athena.hpp"

// M1_INL: forces the body of a KOKKOS_LAMBDA to be inlined into the launch wrapper,
// placed after the parameter list:  KOKKOS_LAMBDA(...) M1_INL {...}.  The large per-cell
// lambdas of the implicit solve are otherwise CALLED from the Kokkos wrapper (the inliner
// refuses them), which pins them at 64 VGPRs and puts a 176-828 byte call frame per lane
// in scratch memory.  Inlining leaves the arithmetic of the lambda as it is.
#define M1_INL __attribute__((always_inline))

namespace radm1 {

template <typename ExeSpace, typename Function>
inline void par_for_lb(const std::string &name, ExeSpace exec_space,
                       const int &nl, const int &nu, const int &kl, const int &ku,
                       const int &jl, const int &ju, const int &il, const int &iu,
                       const Function &function) {
  const int nn = nu - nl + 1;
  const int nk = ku - kl + 1;
  const int nj = ju - jl + 1;
  const int ni = iu - il + 1;
  const int nnkji = nn * nk * nj * ni;
  const int nkji  = nk * nj * ni;
  const int nji   = nj * ni;
  Kokkos::parallel_for(name,
  Kokkos::RangePolicy<ExeSpace, Kokkos::LaunchBounds<256,1>>(exec_space, 0, nnkji),
  KOKKOS_LAMBDA(const int &idx) M1_INL {
    int n = (idx)/nkji;
    int k = (idx - n*nkji)/nji;
    int j = (idx - n*nkji - k*nji)/ni;
    int i = (idx - n*nkji - k*nji - j*ni) + il;
    n += nl;
    k += kl;
    j += jl;
    function(n, k, j, i);
  });
}

} // namespace radm1
#endif // RAD_M1_RAD_M1_PARFOR_HPP_
