//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_memarena.cpp
//! \brief memcut-1113 <rad_m1>/mem_arena: the implicit solve's per-cell work array iw and
//! the vet_gd intensity arrays (VgdRag::d of vgd_i, vgd_itw, vgd_i_xb) share ONE device
//! allocation.  Their lifetimes are disjoint: the intensities are written and read only
//! inside VetGdBuild (sweep, moments, half-range, ring), iw is never read across a
//! VetGdBuild (everything a solve needs from iw is written after the build that precedes
//! it).  Both facts are checked by MEMCUT_POISON (below) and by the bitwise gates.
//!
//! Layout: marena is a (nmb_a, niw, n3, n2, n1) LayoutRight View, nmb_a >= nmb; iw is its
//! first nmb blocks (a contiguous subview: the same type, a managed handle).  The sweep
//! arrays are unmanaged 1-D Views at offsets 0, n_0, n_0 + n_1, ... of the same memory.
//! The allocation grows (copy + rebind) only before the first dt_retry checkpoint
//! (MemArenaSeal from DtRetry::Save), when no other handle of iw or of the intensity
//! arrays exists; afterwards an array that does not fit is allocated on its own, as
//! without the arena.
//!
//! MEMCUT_POISON (diagnostic env, default unset): a comma list of arrays that are filled
//! with NaN (MEMCUT_PVAL=big: +7.7e22) at the entry of every VetGdBuild (iw, ost, thw, thq,
//! krw, odc, fmm_t, vmw, d3w, aa) and at its exit (vgd: every intensity array).  A run whose
//! results stay bitwise with a list poisoned proves the lifetime split for that list.

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "parameter_input.hpp"
#include "rad_m1/rad_m1.hpp"

namespace radm1 {

namespace {
bool EnvListHas(const char *env, const char *name) {
  const char *e = std::getenv(env);
  if (e == nullptr) return false;
  const std::string s = std::string(",") + e + ",";
  return s.find(std::string(",") + name + ",") != std::string::npos ||
         s.find(",all,") != std::string::npos;
}
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::MemArenaIw
//! \brief allocate iw (zero-filled, as Kokkos::realloc did), inside the arena when on

void RadiationM1::MemArenaIw(ParameterInput *pin, const int nmb, const int niw,
                             const int n3, const int n2, const int n1) {
  marena_on = vgd_on;
  if (pin->DoesParameterExist("rad_m1", "mem_arena")) {
    marena_on = marena_on && pin->GetBoolean("rad_m1", "mem_arena");
  }
  const char *e = std::getenv("M1_ARENA");
  if (e != nullptr && std::string(e) == "0") {marena_on = false;}
  if (!marena_on) {
    Kokkos::realloc(iw, nmb, niw, n3, n2, n1);
    return;
  }
  marena = DvceArray5D<Real>("m1_arena", nmb, niw, n3, n2, n1);
  iw = Kokkos::subview(marena, std::make_pair(0, nmb), Kokkos::ALL, Kokkos::ALL,
                       Kokkos::ALL, Kokkos::ALL);
  marena_sw = 0;
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationM1::MemArenaVgd
//! \brief bind a.d (n Reals, zero-filled by the caller) into the arena's sweep part;
//! false: not bound, the caller allocates a.d itself

bool RadiationM1::MemArenaVgd(VgdRag &a, const int64_t n) {
  if (!marena_on || marena.data() == nullptr) return false;
  const int64_t per = static_cast<int64_t>(marena.extent(1))*marena.extent(2)*
                      marena.extent(3)*marena.extent(4);   // Reals per arena block
  const int64_t off = marena_sw;
  const int64_t need = off + n;
  const int64_t have = static_cast<int64_t>(marena.extent(0))*per;
  if (need > have) {
    if (marena_sealed) return false;
    // grow: a new allocation of whole blocks, the old bytes copied in front (iw and the
    // arrays bound so far keep their contents and offsets), then rebind them
    const int nmb_a = static_cast<int>((need + per - 1)/per);
    Kokkos::fence();
    DvceArray5D<Real> nb("m1_arena", nmb_a, marena.extent_int(1), marena.extent_int(2),
                         marena.extent_int(3), marena.extent_int(4));
    Kokkos::deep_copy(Kokkos::subview(nb, std::make_pair(0, marena.extent_int(0)),
                                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL),
                      marena);
    const Real *o0 = marena.data(), *o1 = marena.data() + have;
    Real *n0 = nb.data();
    VgdRag *bound[2] = {&vgd_i, &vgd_itw};   // (+ vgd_i_xb with xffuse-1113)
    for (VgdRag *b : bound) {
      if (b->d.data() != nullptr && b->d.data() >= o0 && b->d.data() < o1) {
        b->d = DvceArray1D<Real>(n0 + (b->d.data() - o0), b->d.extent(0));
      }
    }
    const int nmb = iw.extent_int(0);
    marena = nb;
    iw = Kokkos::subview(marena, std::make_pair(0, nmb), Kokkos::ALL, Kokkos::ALL,
                         Kokkos::ALL, Kokkos::ALL);
  }
  a.d = DvceArray1D<Real>(marena.data() + off, n);
  // 64-Real alignment of the next array
  marena_sw = need + ((64 - (need % 64)) % 64);
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1> mem_arena: intensities " << n*sizeof(Real)/1.0e9
              << " GB in the iw allocation (" << marena.span()*sizeof(Real)/1.0e9
              << " GB, iw " << iw.span()*sizeof(Real)/1.0e9 << " GB) on rank 0"
              << std::endl;
  }
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::MemcutPoison
//! \brief diagnostic (env MEMCUT_POISON, see the file comment): where 0 = VetGdBuild
//! entry (the solve-phase arrays), 1 = its exit (the intensities)

void RadiationM1::MemcutPoison(const int where) {
  static const bool on = (std::getenv("MEMCUT_POISON") != nullptr);
  if (!on) return;
  const char *pv = std::getenv("MEMCUT_PVAL");
  const Real v = (pv != nullptr && std::string(pv) == "big") ? 7.7e22
                 : std::numeric_limits<Real>::quiet_NaN();
  auto p5 = [&](DvceArray5D<Real> &a, const char *nm) {
    if (a.data() != nullptr && EnvListHas("MEMCUT_POISON", nm)) {Kokkos::deep_copy(a, v);}
  };
  if (where == 0) {
    p5(iw, "iw");
    p5(ost, "ost");
    p5(thw, "thw");
    p5(thq, "thq");
    p5(krw, "krw");
    p5(odc, "odc");
    p5(fmm_t, "fmm_t");
    p5(vmw, "vmw");
    p5(d3w, "d3w");
    p5(aa_xc, "aa");
    p5(aa_fc, "aa");
    p5(aa_xp, "aa");
    p5(aa_fp, "aa");
  } else if (EnvListHas("MEMCUT_POISON", "vgd")) {
    VgdRag *bound[2] = {&vgd_i, &vgd_itw};   // (+ vgd_i_xb with xffuse-1113)
    for (VgdRag *b : bound) {
      if (b->d.data() != nullptr) {Kokkos::deep_copy(b->d, v);}
    }
  }
  Kokkos::fence();
}

}  // namespace radm1
