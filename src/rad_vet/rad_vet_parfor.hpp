#ifndef RAD_VET_RAD_VET_PARFOR_HPP_
#define RAD_VET_RAD_VET_PARFOR_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_parfor.hpp
//! \brief 4-D par_for with Kokkos::LaunchBounds<256,2> and flattened 4-D reductions.
//! The launch bound caps registers at 128 on CUDA (occupancy) and lifts the 1024-thread
//! VGPR cap on HIP (no spills); the arithmetic equals athena's par_for.

#include <string>

#include "athena.hpp"

// M1_INL: force-inline a KOKKOS_LAMBDA body into the launch wrapper (no call frame in
// scratch); usage KOKKOS_LAMBDA(...) M1_INL {...}.
#define M1_INL __attribute__((always_inline))

namespace radvet {

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
  Kokkos::RangePolicy<ExeSpace, Kokkos::LaunchBounds<256,2>>(exec_space, 0, nnkji),
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

//! par_reduce_flat2: 4-D (m,k,j,i) parallel_reduce with two reducers in one pass, on a
//! flat 1-D RangePolicy (faster than MDRangePolicy)
template <typename Function, typename Reducer1, typename Reducer2>
inline void par_reduce_flat2(const std::string &name, const int &nl, const int &nu,
                             const int &kl, const int &ku, const int &jl, const int &ju,
                             const int &il, const int &iu, const Function &function,
                             const Reducer1 &reducer1, const Reducer2 &reducer2) {
  const int nn = nu - nl + 1;
  const int nk = ku - kl + 1;
  const int nj = ju - jl + 1;
  const int ni = iu - il + 1;
  const int nnkji = nn * nk * nj * ni;
  const int nkji  = nk * nj * ni;
  const int nji   = nj * ni;
  using VT1 = typename Reducer1::value_type;
  using VT2 = typename Reducer2::value_type;
  Kokkos::parallel_reduce(name, Kokkos::RangePolicy<DevExeSpace>(DevExeSpace(), 0, nnkji),
  KOKKOS_LAMBDA(const int &idx, VT1 &v1, VT2 &v2) {
    int n = (idx)/nkji;
    int k = (idx - n*nkji)/nji;
    int j = (idx - n*nkji - k*nji)/ni;
    int i = (idx - n*nkji - k*nji - j*ni) + il;
    n += nl;
    k += kl;
    j += jl;
    function(n, k, j, i, v1, v2);
  }, reducer1, reducer2);
}

} // namespace radvet
#endif // RAD_VET_RAD_VET_PARFOR_HPP_
