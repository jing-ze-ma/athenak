#ifndef UTILS_SCRATCH_VIEW_HPP_
#define UTILS_SCRATCH_VIEW_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file scratch_view.hpp
//! \brief ScratchView<V>(key, extents...): a persistent scratch View per key, reallocated
//! only when its extents change, so that a per-step temporary costs no cudaMalloc /
//! cudaFree (each cudaFree synchronises the device).  A new allocation is
//! zero-initialised (as a fresh labelled View); a reused one keeps its last contents, so
//! a caller that relies on zeros in elements it does not write zeroes them itself.
//! The Views are released at Kokkos::finalize.

#include <map>
#include <string>

#include <Kokkos_Core.hpp>

template <class V, class... I>
V ScratchView(const std::string &key, const I... ext) {
  static std::map<std::string, V> *m = nullptr;
  if (m == nullptr) {
    m = new std::map<std::string, V>();
    Kokkos::push_finalize_hook([]() {delete m; m = nullptr;});
  }
  V &v = (*m)[key];
  const std::size_t e[] = {static_cast<std::size_t>(ext)...};
  bool same = (v.data() != nullptr);
  for (std::size_t r = 0; r < sizeof...(I) && same; ++r) same = (v.extent(r) == e[r]);
  if (!same) {
    v = V();
    v = V(key, ext...);
  }
  return v;
}

#endif  // UTILS_SCRATCH_VIEW_HPP_
