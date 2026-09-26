#ifndef UTILS_DEEP_COPY_ACROSS_HPP_
#define UTILS_DEEP_COPY_ACROSS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file deep_copy_across.hpp
//! \brief Kokkos::deep_copy between a host View and a device View when either of them is
//! a NON-CONTIGUOUS subview, e.g. one variable n of a (m,n,k,j,i) array over a range of
//! MeshBlocks.
//!
//! Kokkos cannot copy strided data between two memory spaces that no single execution
//! space can reach ("deep_copy with no available copy mechanism"). On a CPU build, or on
//! an APU where the host reaches device memory, a plain deep_copy works, and that is what
//! DeepCopyAcross() does there. Otherwise it stages through contiguous LayoutRight
//! temporaries in each space: strided -> contiguous in the source space, contiguous ->
//! contiguous across, contiguous -> strided in the destination space. A pure copy of the
//! same values either way, so results are unchanged.

#include <string>

#include "athena.hpp"

namespace deep_copy_across {

//! a contiguous LayoutRight View in memory space S with the extents of v
template <class S, class V>
Kokkos::View<typename V::non_const_data_type, Kokkos::LayoutRight, S>
ContigLike(const V &v, const std::string &lbl) {
  constexpr unsigned r = static_cast<unsigned>(V::rank);
  auto e = [&](const unsigned d) -> size_t {
    return (d < r) ? v.extent(d) : KOKKOS_IMPL_CTOR_DEFAULT_ARG;
  };
  Kokkos::LayoutRight lay(e(0), e(1), e(2), e(3), e(4), e(5), e(6), e(7));
  return Kokkos::View<typename V::non_const_data_type, Kokkos::LayoutRight, S>(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, lbl), lay);
}

//! deep_copy(dst, src) that also works for strided views across host and device
template <class DstV, class SrcV>
void DeepCopyAcross(const DstV &dst, const SrcV &src) {
  using DS = typename DstV::memory_space;
  using SS = typename SrcV::memory_space;
  constexpr bool host_reaches_both =
      Kokkos::SpaceAccessibility<Kokkos::HostSpace, DS>::accessible &&
      Kokkos::SpaceAccessibility<Kokkos::HostSpace, SS>::accessible;
  if constexpr (host_reaches_both) {
    Kokkos::deep_copy(dst, src);
  } else {
    auto sc = ContigLike<SS>(src, "dca_src");
    auto dc = ContigLike<DS>(dst, "dca_dst");
    Kokkos::deep_copy(sc, src);
    Kokkos::deep_copy(dc, sc);
    Kokkos::deep_copy(dst, dc);
  }
}

}  // namespace deep_copy_across

#endif  // UTILS_DEEP_COPY_ACROSS_HPP_
