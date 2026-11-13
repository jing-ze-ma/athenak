//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_memarena.cpp
//! \brief memcut-1113 <rad_m1>/mem_arena: the vet_gd intensity arrays (VgdRag::d of
//! vgd_i, vgd_i_xb, vgd_itw) and the implicit solve's scratch arrays that are dead across
//! a VetGdBuild (ost, fmm_t, thw, vmw, thq, odc, d3w, krw, aa_*) share ONE device
//! allocation.  The intensities are written and read only inside VetGdBuild (sweep,
//! moments, half-range, ring); those solve arrays are re-written after every build before
//! they are read.  Both facts are checked by MEMCUT_POISON (below) and the bitwise gates.
//! iw is NOT shared: some of its channels live across a build (poisoned, the run breaks).
//!
//! The intensity arrays are unmanaged 1-D Views at offsets 0, n_0, ... of marena (it
//! grows by copy + rebind while they are allocated, before the first dt_retry
//! checkpoint); at the first checkpoint (MemArenaSeal from DtRetry::Save, before it
//! copies anything) the solve arrays are moved over them (contents copied, members
//! rebound to unmanaged Views at offsets 0, ... of marena), as many as fit.  No
//! rebinding afterwards.
//!
//! MEMCUT_POISON (diagnostic env, default unset): a comma list of arrays that are filled
//! with NaN (MEMCUT_PVAL=big: +7.7e22) at the entry of every VetGdBuild (iw, ost, thw,
//! thq, krw, odc, fmm_t, vmw, d3w, aa) and at its exit (vgd: every intensity array).  A
//! run whose results stay bitwise with a list poisoned shows the lifetime split for it.

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
//! \fn void RadiationM1::MemArenaInit
//! \brief read the key; allocate iw (zero-filled) as before (iw is NOT in the arena: some
//! of its channels live across a VetGdBuild, MEMCUT_POISON=iw, job 13107429)

void RadiationM1::MemArenaInit(ParameterInput *pin, const int nmb, const int niw,
                               const int n3, const int n2, const int n1) {
  marena_on = vgd_on;
  if (pin->DoesParameterExist("rad_m1", "mem_arena")) {
    marena_on = marena_on && pin->GetBoolean("rad_m1", "mem_arena");
  }
  const char *e = std::getenv("M1_ARENA");
  if (e != nullptr && std::string(e) == "0") {marena_on = false;}
  Kokkos::realloc(iw, nmb, niw, n3, n2, n1);
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationM1::MemArenaVgd
//! \brief bind a.d (n Reals, zero-filled by the caller) into the arena (the intensity
//! arrays one after the other); false: not bound, the caller allocates a.d itself

bool RadiationM1::MemArenaVgd(VgdRag &a, const int64_t n) {
  if (!marena_on || marena_sealed) return false;
  const int64_t off = marena_sw;
  const int64_t need = off + n;
  const int64_t have = static_cast<int64_t>(marena.extent(0));
  if (need > have) {
    // grow: a new allocation, the old bytes copied in front (the arrays bound so far keep
    // their contents and offsets), then rebind them (before the first checkpoint only)
    Kokkos::fence();
    DvceArray1D<Real> nb("m1_arena", need);
    if (have > 0) {
      Kokkos::deep_copy(Kokkos::subview(nb, std::make_pair(int64_t(0), have)), marena);
    }
    const Real *o0 = marena.data(), *o1 = marena.data() + have;
    VgdRag *bound[3] = {&vgd_i, &vgd_itw, &vgd_i_xb};
    for (VgdRag *b : bound) {
      if (have > 0 && b->d.data() != nullptr && b->d.data() >= o0 && b->d.data() < o1) {
        b->d = DvceArray1D<Real>(nb.data() + (b->d.data() - o0), b->d.extent(0));
      }
    }
    marena = nb;
  }
  a.d = DvceArray1D<Real>(marena.data() + off, n);
  marena_sw = need + ((64 - (need % 64)) % 64);   // 64-Real alignment of the next one
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::MemArenaSeal
//! \brief called by every DtRetry::Save before it copies; the first call moves the solve
//! scratch arrays that are dead across a VetGdBuild (MEMCUT_POISON: ost, fmm_t, thw, vmw,
//! thq, odc, d3w, krw, aa_*, job 13105648) into the arena, over the intensities (contents
//! copied, member rebound), as many as fit; then no more rebinding

void RadiationM1::MemArenaSeal() {
  if (marena_sealed) return;
  marena_sealed = true;
  if (!marena_on || marena.data() == nullptr) return;
  DvceArray5D<Real> *cand[] = {&ost, &fmm_t, &thw, &vmw, &thq, &odc, &d3w, &krw,
                               &aa_xc, &aa_fc, &aa_xp, &aa_fp};
  const char *nm[] = {"ost", "fmm_t", "thw", "vmw", "thq", "odc", "d3w", "krw",
                      "aa_xc", "aa_fc", "aa_xp", "aa_fp"};
  const int64_t have = static_cast<int64_t>(marena.extent(0));
  int64_t off = 0;
  std::string moved;
  Kokkos::fence();
  for (int q = 0; q < static_cast<int>(sizeof(cand)/sizeof(cand[0])); ++q) {
    DvceArray5D<Real> &a = *cand[q];
    const int64_t n = static_cast<int64_t>(a.span());
    if (a.data() == nullptr || n == 0 || !a.span_is_contiguous() || off + n > have) {
      continue;
    }
    DvceArray5D<Real> b(marena.data() + off, a.extent(0), a.extent(1), a.extent(2),
                        a.extent(3), a.extent(4));
    Kokkos::deep_copy(b, a);
    a = b;
    off += n + ((64 - (n % 64)) % 64);
    moved += std::string(" ") + nm[q];
  }
  Kokkos::fence();
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1> mem_arena: intensities " << marena_sw*sizeof(Real)/1.0e9
              << " GB in one allocation; solve scratch moved over them (" << off*
              sizeof(Real)/1.0e9 << " GB on rank 0):" << moved << std::endl;
  }
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
    VgdRag *bound[3] = {&vgd_i, &vgd_itw, &vgd_i_xb};
    for (VgdRag *b : bound) {
      if (b->d.data() != nullptr) {Kokkos::deep_copy(b->d, v);}
    }
  }
  Kokkos::fence();
}

}  // namespace radm1
