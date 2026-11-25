//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file restart.cpp
//! \brief writes restart files

#include <sys/stat.h>  // mkdir

#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <cstdio>      // fwrite(), fclose(), fopen(), fnprintf(), snprintf()
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <type_traits>
#include <utility> // make_pair
#include <vector>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "geodesic-grid/geodesic_grid.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "hydro/hydro_am.hpp"
#include "mhd/mhd.hpp"
#include "coordinates/adm.hpp"
#include "z4c/compact_object_tracker.hpp"
#include "z4c/z4c.hpp"
#include "radiation/radiation.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "srcterms/turb_driver.hpp"
#include "pgen/pgen.hpp"
#include "utils/two_stream_warm_rst.hpp"
#include "utils/two_stream_ck_rst.hpp"
#include "utils/deep_copy_across.hpp"

using deep_copy_across::DeepCopyAcross;
//#include "outputs.hpp"

//----------------------------------------------------------------------------------------
// constructor: also calls BaseTypeOutput base class constructor

RestartOutput::RestartOutput(ParameterInput *pin, Mesh *pm, OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  // create directories for outputs. Comments in binary.cpp constructor explain why
  win_mb_ = pin->GetOrAddReal("job", "restart_window_mb", 1024.0);
  mkdir("rst",0775);
  bool single_file_per_rank = op.single_file_per_rank;
  if (single_file_per_rank) {
    char rank_dir[20];
    std::snprintf(rank_dir, sizeof(rank_dir), "rst/rank_%08d/", global_variable::my_rank);
    mkdir(rank_dir, 0775);
  }
}

//----------------------------------------------------------------------------------------
// RestartOutput::LoadOutputData()
// overload of standard load data function specific to restarts.  Loads dependent
// variables, including ghost zones.

void RestartOutput::LoadOutputData(Mesh *pm) {
  // get spatial dimensions of arrays, including ghost zones
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int nout1 = indcs.nx1 + 2*(indcs.ng);
  int nout2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int nout3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  int nmb = pm->pmb_pack->nmb_thispack;

  // calculate total number of CC variables
  hydro::Hydro* phydro = pm->pmb_pack->phydro;
  mhd::MHD* pmhd = pm->pmb_pack->pmhd;
  adm::ADM* padm = pm->pmb_pack->padm;
  z4c::Z4c* pz4c = pm->pmb_pack->pz4c;
  radiation::Radiation* prad = pm->pmb_pack->prad;
  // <rad_vet>: the four evolved moments (E, F_i), written and read as one more
  // cell-centred block.  Guarded by the pointer everywhere, so a file written by a run
  // WITHOUT <rad_vet> has exactly the layout and the bytes it had before (design sect.6)
  radvet::RadiationVet* pradvet = pm->pmb_pack->pradvet;
  TurbulenceDriver* pturb=pm->pmb_pack->pturb;
  int nhydro=0, nmhd=0, nrad=0, nm1=0, nforce=3, nadm=0, nz4c=0;
  if (pradvet != nullptr) {
    nm1 = radvet::VET_NVAR;
  }
  if (phydro != nullptr) {
    nhydro = phydro->nhydro + phydro->nscalars;
  }
  if (pmhd != nullptr) {
    nmhd = pmhd->nmhd + pmhd->nscalars;
  }
  if (pz4c != nullptr) {
    nz4c = pz4c->nz4c;
  } else if (padm != nullptr) {
    nadm = padm->nadm;
  }
  // if the spacetime is evolved, we do not need to checkpoint/recover the ADM variables
  if (prad != nullptr) {
    nrad = prad->prgeo->nangles;
  }

  // Note for restarts, outarrays are dimensioned (m,n,k,j,i)
  // where the host reads device memory (MI300A, CPU builds) and the dump is written
  // synchronously, a staging array is the contiguous device array (or its device pack)
  // itself instead of a host copy of it: the same bytes, one copy fewer at the dump
  constexpr bool um_ = Kokkos::SpaceAccessibility<HostMemSpace, DevMemSpace>::accessible;
  const bool zc = um_;
  rst_keep_.clear();
  auto stage = [&](auto &out, const auto &sv) {
    using H = std::remove_reference_t<decltype(out)>;
    if (zc && sv.span_is_contiguous()) {
      if constexpr (H::rank == 5) {
        out = H(sv.data(), sv.extent(0), sv.extent(1), sv.extent(2), sv.extent(3),
                sv.extent(4));
      } else {
        out = H(sv.data(), sv.extent(0), sv.extent(1), sv.extent(2), sv.extent(3));
      }
    } else {
      if constexpr (H::rank == 5) {
        Kokkos::realloc(out, sv.extent(0), sv.extent(1), sv.extent(2), sv.extent(3),
                        sv.extent(4));
      } else {
        Kokkos::realloc(out, sv.extent(0), sv.extent(1), sv.extent(2), sv.extent(3));
      }
      Kokkos::deep_copy(out, sv);
    }
  };
  auto hold = [&](HostArray5D<Real> &out, const DvceArray5D<Real> &tmp) {
    if (zc) rst_keep_.push_back(tmp);
    stage(out, tmp);
  };
  // hgmem: a module's pack (fn writes a device-typed (n0, nf, k, j, i) View) where the
  // host does not read device memory: straight into pinned host memory (no device copy
  // at a dump; the same bytes), then staged; elsewhere the device array as before
  auto pack = [&](HostArray5D<Real> &out, const char *lbl, const int n0, const int nf,
                  const auto &fn) {
    if constexpr (!um_) {
      Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> pin(std::string(lbl),
          static_cast<size_t>(n0)*nf*nout3*nout2*nout1);
      DvceArray5D<Real> d(pin.data(), n0, nf, nout3, nout2, nout1);
      fn(d);
      Kokkos::fence();
      stage(out, HostArray5D<Real>(pin.data(), n0, nf, nout3, nout2, nout1));
    } else {
      DvceArray5D<Real> tmp(std::string(lbl), n0, nf, nout3, nout2, nout1);
      fn(tmp);
      hold(out, tmp);
    }
  };
  if (phydro != nullptr) {
    stage(outarray_hyd, Kokkos::subview(phydro->u0, std::make_pair(0,nmb), Kokkos::ALL,
                                        Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (pmhd != nullptr) {
    Kokkos::realloc(outarray_mhd, nmb, nmhd, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_mhd, Kokkos::subview(pmhd->u0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    Kokkos::realloc(outfield.x1f, nmb, nout3, nout2, nout1+1);
    Kokkos::deep_copy(outfield.x1f, Kokkos::subview(pmhd->b0.x1f, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    Kokkos::realloc(outfield.x2f, nmb, nout3, nout2+1, nout1);
    Kokkos::deep_copy(outfield.x2f, Kokkos::subview(pmhd->b0.x2f, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    Kokkos::realloc(outfield.x3f, nmb, nout3+1, nout2, nout1);
    Kokkos::deep_copy(outfield.x3f, Kokkos::subview(pmhd->b0.x3f, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  // the general-EOS temperature cache, see the note on outarray_wth in outputs.hpp
  if (phydro != nullptr && phydro->peos->eos_data.IsGeneral()) {
    stage(outarray_wth, Kokkos::subview(phydro->wtemp, std::make_pair(0,nmb), Kokkos::ALL,
                                        Kokkos::ALL, Kokkos::ALL));
  }
  if (pmhd != nullptr && pmhd->peos->eos_data.IsGeneral()) {
    Kokkos::realloc(outarray_wtm, nmb, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_wtm, Kokkos::subview(pmhd->wtemp,
                      std::make_pair(0,nmb), Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  // the general-EOS derived cache p and Gamma_1, see the note on outarray_wdp
  if (phydro != nullptr && phydro->peos->eos_data.IsGeneral()) {
    Kokkos::realloc(outarray_wdp, nmb, nout3, nout2, nout1);
    DeepCopyAcross(outarray_wdp, Kokkos::subview(phydro->wder,
                   std::make_pair(0,nmb), static_cast<int>(IDPR),
                   Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    Kokkos::realloc(outarray_wdg, nmb, nout3, nout2, nout1);
    DeepCopyAcross(outarray_wdg, Kokkos::subview(phydro->wder,
                   std::make_pair(0,nmb), static_cast<int>(IDG1),
                   Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  // and the MHD one, which is the same cache of the same two channels: see the note on
  // outarray_wdpm.  The general-EOS MHD runs (the deep hot Jupiter) are chained restarts
  // and every link boundary was a discontinuity without it.
  if (pmhd != nullptr && pmhd->peos->eos_data.IsGeneral()) {
    Kokkos::realloc(outarray_wdpm, nmb, nout3, nout2, nout1);
    DeepCopyAcross(outarray_wdpm, Kokkos::subview(pmhd->wder,
                   std::make_pair(0,nmb), static_cast<int>(IDPR),
                   Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    Kokkos::realloc(outarray_wdgm, nmb, nout3, nout2, nout1);
    DeepCopyAcross(outarray_wdgm, Kokkos::subview(pmhd->wder,
                   std::make_pair(0,nmb), static_cast<int>(IDG1),
                   Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  // the internal energy the last inversion produced, see kEintRstMagic in pgen.hpp
  if (phydro != nullptr && phydro->RstEint()) {   // + ideal gas
    stage(outarray_wpv, Kokkos::subview(phydro->w0, std::make_pair(0,nmb), Kokkos::ALL,
                                        Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    Kokkos::realloc(outarray_weh, nmb, nout3, nout2, nout1);
    DeepCopyAcross(outarray_weh, Kokkos::subview(phydro->w0,
                   std::make_pair(0,nmb), static_cast<int>(IEN),
                   Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (pmhd != nullptr && pmhd->RstEint()) {   // + ideal gas
    Kokkos::realloc(outarray_wem, nmb, nout3, nout2, nout1);
    DeepCopyAcross(outarray_wem, Kokkos::subview(pmhd->w0,
                   std::make_pair(0,nmb), static_cast<int>(IEN),
                   Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  // the implicit correlated-k solver's cross-call state (utils/two_stream_ck_rst.hpp)
  ck_rst_hdr_bytes.clear();
  if (two_stream_rt::ck_rst_collect_fn != nullptr) {
    two_stream_rt::CkRstHdr ckh;
    const int nck = (two_stream_rt::ck_rst_collect_fn)(pm, outarray_ck, ckh, nmb, nout3,
                                                       nout2, nout1);
    if (nck > 0) {
      ck_rst_hdr_bytes.resize(sizeof(ckh));
      std::memcpy(ck_rst_hdr_bytes.data(), &ckh, sizeof(ckh));
    }
  }
  // the mode-3 Newton warm-start history, see the note on outarray_wm1
  {
    const int nwm = two_stream_rt::RtWarmLevels();
    if (nwm > 0) {
      Kokkos::realloc(outarray_wm1, nmb, nout3, nout2, nout1);
      Kokkos::deep_copy(outarray_wm1, Kokkos::subview(*two_stream_rt::rt_c3bp_ptr,
                        std::make_pair(0,nmb), Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    }
    if (nwm > 1) {
      Kokkos::realloc(outarray_wm2, nmb, nout3, nout2, nout1);
      Kokkos::deep_copy(outarray_wm2, Kokkos::subview(*two_stream_rt::rt_c3bp2_ptr,
                        std::make_pair(0,nmb), Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    }
  }
  if (prad != nullptr) {
    Kokkos::realloc(outarray_rad, nmb, nrad, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_rad, Kokkos::subview(prad->i0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (pradvet != nullptr) {
    stage(outarray_m1, Kokkos::subview(pradvet->u0, std::make_pair(0,nmb), Kokkos::ALL,
                                       Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    // milestone 3a: the persistent x1 face fluxes of the implicit scheme
    {
      stage(outarray_m1f, Kokkos::subview(pradvet->f0x1, std::make_pair(0,nmb),
                                          Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    }
    // milestone 3b phase B: the transverse face fluxes, in the SAME order the reader
    // (src/pgen/pgen.cpp) takes them back: f0x1, then f0x2, then f0x3.
    if (pradvet->trans_on) {
      stage(outarray_m1f2, Kokkos::subview(pradvet->f0x2, std::make_pair(0,nmb),
                                           Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
      if (pradvet->trans_x3) {
        stage(outarray_m1f3, Kokkos::subview(pradvet->f0x3, std::make_pair(0,nmb),
                                             Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
      }
    }
    // implicit_predictor = step: the stored previous-step increment (radvet::ipred)
    if (pradvet->impl_pred && pradvet->pred_ok) {
      const int npc = static_cast<int>(pradvet->ipred.extent(1));
      Kokkos::realloc(outarray_m1p, nmb, npc, nout3, nout2, nout1);
      Kokkos::deep_copy(outarray_m1p, Kokkos::subview(pradvet->ipred,
                        std::make_pair(0,nmb), Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                        Kokkos::ALL));
    }
    // time_scheme = hesdirk2: the FSAL slope and what goes with it
    if (pradvet->Time2RstNch() > 0) {
      const int nct = pradvet->Time2RstNch();
      if (pradvet->rst_pack_mb) {
        // restart_pack_per_mb: one MeshBlock at a time (same bytes); the clean rays'
        // intensities (the last nfv channels) are copied to the host directly
        const int nfv = pradvet->cln_cssc ? 0 : pradvet->CleanRstNfv(), nc0 = nct - nfv;
        // (the staging of a MeshBlock in pinned host memory, which the device writes
        // directly: no device copy of its channels at a dump; where the host reads device
        // memory (MI300A, CPU builds) the device array itself, packed in place: no
        // staging, no host copy, no host memset)
        constexpr bool um =
            Kokkos::SpaceAccessibility<HostMemSpace, DevMemSpace>::accessible;
        // (independent I/O, synchronous: one MeshBlock's buffer, packed and written
        // MeshBlock by MeshBlock in WriteOutputFile)
        const int nst = (um && nfv == 0 && !IOWrapper::collective) ? 1 : nmb;
        // (and, without ray intensities in the record, in windows of at most
        // <job>/restart_window_mb)
        const int nwin = (nst < nmb && pradvet->CleanRstNfv() == 0) ?
            std::max(1, std::min(nct, static_cast<int>(win_mb_*1048576.0/
                                 (sizeof(Real)*nout1*nout2*nout3)))) : nct;
        if (um) {
          outarray_m1t_d = DvceArray5D<Real>("rst-m1t-d", nst, nwin, nout3, nout2, nout1);
          outarray_m1t = HostArray5D<Real>(outarray_m1t_d.data(), nst, nwin, nout3, nout2,
                                           nout1);
        } else {
          Kokkos::realloc(outarray_m1t, nmb, nct, nout3, nout2, nout1);
        }
        Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> pin("rst-m1t-mb",
            um ? 0 : static_cast<size_t>(nc0)*nout3*nout2*nout1);
        auto AL = Kokkos::ALL;
        for (int m = 0; m < ((nst < nmb) ? 0 : nmb); ++m) {
          DvceArray5D<Real> tmp(um ? outarray_m1t_d.data() + outarray_m1t_d.stride(0)*m :
                                pin.data(), 1, nc0, nout3, nout2, nout1);
          pradvet->Time2RstPack(tmp, 1, m);
          if (!um) {
            Kokkos::deep_copy(Kokkos::subview(outarray_m1t, m, std::make_pair(0, nc0), AL,
                                              AL, AL),
                              Kokkos::subview(tmp, 0, AL, AL, AL, AL));
          }
          auto fv = Kokkos::subview(outarray_m1t, m, std::make_pair(nc0, nct), AL, AL,
                                    AL);
          if (nfv > 0 && pradvet->CleanRstFv().extent_int(0) > m) {
            Kokkos::deep_copy(fv, Kokkos::subview(pradvet->CleanRstFv(), m, AL, AL, AL,
                                                  AL));
          } else if (nfv > 0) {
            Kokkos::deep_copy(fv, 0.0);
          }
        }
        if (um) Kokkos::fence();
      } else {
        DvceArray5D<Real> tmp("rst-m1t", nmb, nct, nout3, nout2, nout1);
        pradvet->Time2RstPack(tmp, nmb);
        Kokkos::realloc(outarray_m1t, nmb, nct, nout3, nout2, nout1);
        Kokkos::deep_copy(outarray_m1t, tmp);
      }
    }
    // implicit_closure_thin_relax: the relaxed closure of the last step (ctr_mem)
    if (pradvet->impl_ctrelax > 0.0 && pradvet->ctr_init) {
      Kokkos::realloc(outarray_m1c, nmb, 4, nout3, nout2, nout1);
      Kokkos::deep_copy(outarray_m1c, Kokkos::subview(pradvet->ctr_mem,
                        std::make_pair(0,nmb), Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                        Kokkos::ALL));
    }
    // <hydro|mhd>/rad_signal_speed: the signal-speed inputs of the last M1 step
    if ((phydro != nullptr && phydro->rad_signal_speed) ||
        (pmhd != nullptr && pmhd->rad_signal_speed)) {
      const int rmode = pradvet->RssRstMode();
      const int nrs = pradvet->RssRstNch(rmode);
      // (zero-copy, independent I/O: one MeshBlock's pack, filled MeshBlock by MeshBlock
      // by the writer)
      const bool one = zc && nmb > 1 && !IOWrapper::collective;
      pack(outarray_m1s, "rst-m1s", one ? 1 : nmb, nrs, [&](DvceArray5D<Real> &t) {
        if (!one) pradvet->RssRstPack(t, nmb, rmode);
      });
    }
    // implicit_gas_src = conservative(2): the cell-F slope csk1
    if (pradvet->CskRstNch() > 0) {
      Kokkos::realloc(outarray_m1k, nmb, 3, nout3, nout2, nout1);
      Kokkos::deep_copy(outarray_m1k, Kokkos::subview(pradvet->csk1,
                                                      std::make_pair(0,nmb),
                        Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    }
    // vet_rotate_avg: the averaging history, densely packed slabs
    const int nra = pradvet->RotAvgRstNslab(nout1*nout2*nout3);
    if (nra > 0) {
      // packed per MeshBlock in chunks of 4 slabs (was a device copy of the
      // whole record, the largest transient on top of the steady state); same bytes
      Kokkos::realloc(outarray_m1r, nmb, nra, nout3, nout2, nout1);
      pradvet->RotAvgRstPackHost(outarray_m1r, nmb, 4);
    }
  }
  // <hydro|mhd>/am_vector: the spin reservoir (+ the radiation's hesdirk2
  // spin slope amk1), kAmSpinRstMagic, the last slabs of each record
  {
    hydro::AmModule *pam = hydro::AmModuleOf(pm->pmb_pack);
    const bool k1 = (pam != nullptr) && pam->am_vector && (pradvet != nullptr) &&
                    (pradvet->amk1.extent_int(0) > 0);
    const int nams = ((pam != nullptr) && pam->am_vector) ? (k1 ? 2 : 1) : 0;
    if (nams > 0) {
      auto mb = std::make_pair(0, nmb);
      auto AL = Kokkos::ALL;
      DvceArray5D<Real> tmp("rst-ams", nmb, nams, nout3, nout2, nout1);
      Kokkos::deep_copy(Kokkos::subview(tmp, AL, 0, AL, AL, AL),
                        Kokkos::subview(pam->am_spin, mb, AL, AL, AL));
      if (k1) {
        Kokkos::deep_copy(Kokkos::subview(tmp, AL, 1, AL, AL, AL),
                          Kokkos::subview(pradvet->amk1, mb, AL, AL, AL));
      }
      hold(outarray_ams, tmp);
    }
  }
  // the per-cell pgen slabs (kPgenSlabRstMagic), behind the spin ones
  {
    DvceArray5D<Real> *ps = (pm->pgen != nullptr) ? pm->pgen->pgen_rst_slabs : nullptr;
    const int nps = (ps != nullptr) ? ps->extent_int(1) : 0;
    if (nps > 0) {
      DvceArray5D<Real> tmp("rst-pgs", nmb, nps, nout3, nout2, nout1);
      Kokkos::deep_copy(tmp, Kokkos::subview(*ps, std::make_pair(0, nmb), Kokkos::ALL,
                                             Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
      hold(outarray_pgs, tmp);
    }
  }
  // implicit_hr_cmem_rst: the applied face fluxes (kVetFaxRstMagic)
  if (pradvet != nullptr && pradvet->FaxRstNslab() > 0) {
    const int nf = pradvet->FaxRstNslab();
    pack(outarray_fax, "rst-fax", nmb, nf, [&](DvceArray5D<Real> &t) {
      pradvet->FaxRstPack(t, nmb);
    });
  }
  // vet_x1_periodic: the warm-start exit plane (kVetX1eRstMagic)
  if (pradvet != nullptr && pradvet->X1eRstNslab(nout1) > 0) {
    pack(outarray_x1e, "rst-x1e", nmb, pradvet->X1eRstNslab(nout1),
         [&](DvceArray5D<Real> &t) {pradvet->X1eRstPack(t, nmb);});
  }
  if (pturb != nullptr) {
    Kokkos::realloc(outarray_force, nmb, nforce, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_force, Kokkos::subview(pturb->force, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (pz4c != nullptr) {
    Kokkos::realloc(outarray_z4c, nmb, nz4c, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_z4c, Kokkos::subview(pz4c->u0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  } else if (padm != nullptr) {
    Kokkos::realloc(outarray_adm, nmb, nadm, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_adm, Kokkos::subview(padm->u_adm, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }

  // calculate max/min number of MeshBlocks across all ranks
  noutmbs_max = pm->nmb_eachrank[0];
  noutmbs_min = pm->nmb_eachrank[0];
  for (int i=0; i<(global_variable::nranks); ++i) {
    noutmbs_max = std::max(noutmbs_max,pm->nmb_eachrank[i]);
    noutmbs_min = std::min(noutmbs_min,pm->nmb_eachrank[i]);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RestartOutput:::WriteOutputFile(Mesh *pm)
//  \brief Cycles over all MeshBlocks and writes everything to a single restart file

void RestartOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  // get spatial dimensions of arrays, including ghost zones
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int nout1 = indcs.nx1 + 2*(indcs.ng);
  int nout2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int nout3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  hydro::Hydro* phydro = pm->pmb_pack->phydro;
  mhd::MHD* pmhd = pm->pmb_pack->pmhd;
  radiation::Radiation* prad = pm->pmb_pack->prad;
  radvet::RadiationVet* pradvet = pm->pmb_pack->pradvet;
  TurbulenceDriver* pturb=pm->pmb_pack->pturb;
  z4c::Z4c* pz4c = pm->pmb_pack->pz4c;
  adm::ADM* padm = pm->pmb_pack->padm;
  // the OPTIONAL, MARKED pgen state block (see ProblemGenerator::pgen_rst_write_func).
  // Empty -- which is every problem generator but red_giant with problem/mlt_alpha > 0 --
  // and not one byte is written, so the file layout is exactly what it always was.
  std::vector<char> pgen_state;
  if (pm->pgen != nullptr && pm->pgen->pgen_rst_write_func != nullptr) {
    pgen_state = (pm->pgen->pgen_rst_write_func)();
  }
  // the mode-3 Newton warm-start history: how many levels this run keeps (0 when the
  // warm start is off, or on but not yet allocated), and the little marked header that
  // says so -- int32 nlev, int32 pad, then the previous call's bdt, which is the
  // denominator of the rt_impl_warm = 2 extrapolation ratio and is state as much as the
  // levels are.  See utils/two_stream_warm_rst.hpp.
  const int nwarm = two_stream_rt::RtWarmLevels();
  char warm_hdr[2*sizeof(std::int32_t) + sizeof(Real)];
  {
    const std::int32_t hdr[2] = {static_cast<std::int32_t>(nwarm), 0};
    const Real bdtp = two_stream_rt::rt_c3_bdt_prev;
    std::memcpy(&(warm_hdr[0]), &(hdr[0]), sizeof(hdr));
    std::memcpy(&(warm_hdr[0]) + sizeof(hdr), &bdtp, sizeof(bdtp));
  }
  // <rad_vet> implicit_predictor = step: the number of ipred channels written (0 when no
  // increment is stored, and then not one byte), and its marked header -- int32 have,
  // int32 ncomp, Real pred_dt.  See radvet::kVetPredRstMagic.
  const int npred = (pradvet != nullptr && pradvet->impl_pred && pradvet->pred_ok) ?
                    static_cast<int>(pradvet->ipred.extent(1)) : 0;
  char pred_hdr[2*sizeof(std::int32_t) + sizeof(Real)];
  {
    const std::int32_t hdr[2] = {(npred > 0) ? 1 : 0, static_cast<std::int32_t>(npred)};
    const Real pdt = (pradvet != nullptr) ? pradvet->pred_dt : 0.0;
    std::memcpy(&(pred_hdr[0]), &(hdr[0]), sizeof(hdr));
    std::memcpy(&(pred_hdr[0]) + sizeof(hdr), &pdt, sizeof(pdt));
  }
  // <rad_vet> time_scheme = hesdirk2: the slope block, behind the predictor one --
  // int32 have, int32 nch, Real pred2_dt, Real dt_prev, Real vprev
  const int nt2 = (pradvet != nullptr) ? pradvet->Time2RstNch() : 0;
  // <rad_vet> implicit_one_pass: its state, behind the hesdirk2 header
  // (9 Reals; 18 with implicit_one_pass_auto, whose switch-off state follows)
  const bool wonep = (pradvet != nullptr) && (pradvet->impl_onep > 0);
  Real onep_hdr[18];
  IOWrapperSizeT onep_nb = 0;
  if (wonep) {
    for (int t = 0; t < 3; ++t) {
      onep_hdr[t] = pradvet->onep_qa[t];
      onep_hdr[3+t] = pradvet->onep_qb[t];
      onep_hdr[6+t] = pradvet->onep_cnt[t];
      onep_hdr[9+t] = pradvet->onep_off[t];
      onep_hdr[12+t] = pradvet->onep_actr[t];
      onep_hdr[15+t] = pradvet->onep_prb[t];
    }
    onep_nb = (pradvet->impl_onep_auto ? 18 : 9)*sizeof(Real);
  }
  // the time2_auto_be back-off state, only while a skip run is active
  const bool wabe = (pradvet != nullptr) && (pradvet->t2_abe > 0) &&
                    (pradvet->t2_abe_len > 0 || pradvet->t2_abe_left > 0);
  std::int32_t abe_hdr[2] = {0, 0};
  if (wabe) {
    abe_hdr[0] = static_cast<std::int32_t>(pradvet->t2_abe_len);
    abe_hdr[1] = static_cast<std::int32_t>(pradvet->t2_abe_left);
  }
  char t2_hdr[2*sizeof(std::int32_t) + 3*sizeof(Real)];
  {
    // (version 2: the staggered fkin block carries K1_d; version 3, the
    // clean scheme: hv[2] also holds t2_ok and cln_k1ok, as bits 2 and 4; versions 4
    // and 5, the clean scheme: no zero channels, the ring dense, then also the cubed
    // sphere's sc rays, Time2RstNchWant; version 7: the ring without the cubed sphere's
    // zero stress / CLN_G fields, CleanRingFields)
    const bool cl = (pradvet != nullptr);
    const bool fk = cl && pradvet->cln_fsolve && pradvet->cln_fkin;
    const std::int32_t hdr[2] = {(nt2 > 0) ? (cl ? 7 : (fk ? 2 : 1)) : 0,
                                 static_cast<std::int32_t>(nt2)};
    Real hv[3] = {0.0, 0.0, 0.0};
    if (pradvet != nullptr) {
      hv[0] = pradvet->pred2_ok ? pradvet->pred2_dt : 0.0;
      hv[1] = pradvet->t2_dtprev;
      hv[2] = (pradvet->t2_vprev ? 1.0 : 0.0) + (cl && pradvet->t2_ok ? 2.0 : 0.0) +
              (cl && pradvet->cln_k1ok ? 4.0 : 0.0);
    }
    std::memcpy(&(t2_hdr[0]), &(hdr[0]), sizeof(hdr));
    std::memcpy(&(t2_hdr[0]) + sizeof(hdr), &(hv[0]), sizeof(hv));
  }
  // the general-EOS internal energy slabs (kEintRstMagic): one per general-EOS module
  // and every ideal-gas module (RstEint, <hydro|mhd>/rst_eint)
  const bool we_hyd = (phydro != nullptr) && phydro->RstEint();
  const bool we_mhd = (pmhd != nullptr) && pmhd->RstEint();
  const int neint = (we_hyd ? 1 : 0) + (we_mhd ? 1 : 0);
  // the implicit-ck state block (utils/two_stream_ck_rst.hpp), behind even that one
  const int nck = ck_rst_hdr_bytes.empty() ? 0
                  : static_cast<int>(outarray_ck.extent(0));
  // <rad_vet> implicit_closure_thin_relax: 4 ctr_mem slabs behind everything (0 = none,
  // and then not one byte).  radvet::kVetCtrRstMagic.
  const int nctr = (pradvet != nullptr && pradvet->impl_ctrelax > 0.0 &&
                    pradvet->ctr_init)
                   ? 4 : 0;
  char ctr_hdr[2*sizeof(std::int32_t)];
  {
    const std::int32_t hdr[2] = {(nctr > 0) ? 1 : 0, static_cast<std::int32_t>(nctr)};
    std::memcpy(&(ctr_hdr[0]), &(hdr[0]), sizeof(hdr));
  }
  // <hydro|mhd>/rad_signal_speed: nrss slabs behind even those (radvet::kVetRssRstMagic)
  const bool wrss = (pradvet != nullptr) &&
                    ((phydro != nullptr && phydro->rad_signal_speed) ||
                     (pmhd != nullptr && pmhd->rad_signal_speed));
  const int rss_mode = wrss ? pradvet->RssRstMode() : 0;
  const int nrss = wrss ? pradvet->RssRstNch(rss_mode) : 0;
  char rss_hdr[4*sizeof(std::int32_t)];
  {
    const std::int32_t hdr[4] = {(nrss > 0) ? 1 : 0, static_cast<std::int32_t>(nrss),
                                 static_cast<std::int32_t>(rss_mode), 0};
    std::memcpy(&(rss_hdr[0]), &(hdr[0]), sizeof(hdr));
  }
  // <rad_vet>/implicit_gas_src = conservative(2): csk1 + the flags,
  // ncsk slabs behind the vet_gd_async ones (radvet::kVetCskRstMagic)
  const bool wcsk = (pradvet != nullptr) && pradvet->CskRstHdr();
  const int ncsk = wcsk ? pradvet->CskRstNch() : 0;
  char csk_hdr[4*sizeof(std::int32_t)];
  {
    const std::int32_t hdr[4] = {wcsk ? 1 : 0, static_cast<std::int32_t>(ncsk),
                                 (wcsk && pradvet->csk1_ok) ? 1 : 0, 0};
    std::memcpy(&(csk_hdr[0]), &(hdr[0]), sizeof(hdr));
  }
  // <rad_vet>/vet_rotate_avg: the averaging history, nra slabs behind
  // even the csk ones (radvet::kVetRaRstMagic; RadiationVet::RotAvgRstHdr)
  std::int64_t ra_hdr[radvet::RadiationVet::VET_RA_RST_NH] = {0};
  if (pradvet != nullptr) {pradvet->RotAvgRstHdr(&(ra_hdr[0]));}
  const int nra = static_cast<int>(ra_hdr[1]);
  // am_vector spin slabs (hydro::kAmSpinRstMagic), behind even the nra ones
  int nams = 0;
  {
    hydro::AmModule *pam = hydro::AmModuleOf(pm->pmb_pack);
    if ((pam != nullptr) && pam->am_vector) {
      nams = ((pradvet != nullptr) && (pradvet->amk1.extent_int(0) > 0)) ? 2 : 1;
    }
  }
  char ams_hdr[4*sizeof(std::int32_t)];
  {
    const std::int32_t hdr[4] = {1, static_cast<std::int32_t>(nams), 0, 0};
    std::memcpy(&(ams_hdr[0]), &(hdr[0]), sizeof(hdr));
  }
  // the per-cell pgen slabs, behind even the spin ones
  const int npgs = (pm->pgen != nullptr && pm->pgen->pgen_rst_slabs != nullptr) ?
                   pm->pgen->pgen_rst_slabs->extent_int(1) : 0;
  char pgs_hdr[4*sizeof(std::int32_t)];
  {
    const std::int32_t hdr[4] = {1, static_cast<std::int32_t>(npgs), 0, 0};
    std::memcpy(&(pgs_hdr[0]), &(hdr[0]), sizeof(hdr));
  }
  // the rad_vet applied face fluxes, behind even the pgen slabs
  const int nfax = (pradvet != nullptr) ? pradvet->FaxRstNslab() : 0;
  char fax_hdr[4*sizeof(std::int32_t)];
  {
    const std::int32_t hdr[4] = {1, static_cast<std::int32_t>(nfax), 0, 0};
    std::memcpy(&(fax_hdr[0]), &(hdr[0]), sizeof(hdr));
  }
  // the vet_x1_periodic warm-start plane, behind the fax slabs
  const int nx1e = (pradvet != nullptr) ? pradvet->X1eRstNslab(nout1) : 0;
  const std::int32_t x1e_hdr[4] = {1, static_cast<std::int32_t>(nx1e), 0, 0};
  char eint_hdr[2*sizeof(std::int32_t)];
  {
    // (hdr[1]: the hydro w0 density and velocity slabs behind the hydro w0(IEN))
    const std::int32_t hdr[2] = {static_cast<std::int32_t>(neint), we_hyd ? 4 : 0};
    std::memcpy(&(eint_hdr[0]), &(hdr[0]), sizeof(hdr));
  }
  int nhydro=0, nmhd=0, nrad=0, nm1=0, nforce=3, nz4c=0, nadm=0, nco=0;
  if (pradvet != nullptr) {
    nm1 = radvet::VET_NVAR;
  }
  if (phydro != nullptr) {
    nhydro = phydro->nhydro + phydro->nscalars;
  }
  if (pmhd != nullptr) {
    nmhd = pmhd->nmhd + pmhd->nscalars;
  }
  if (prad != nullptr) {
    nrad = prad->prgeo->nangles;
  }
  if (pz4c != nullptr) {
    nz4c = pz4c->nz4c;
    nco = pz4c->ptracker.size();
  } else if (padm != nullptr) {
    nadm = padm->nadm;
  }
  bool single_file_per_rank = out_params.single_file_per_rank;
  std::string fname;
  if (single_file_per_rank) {
    // Generate a directory and filename for each rank
    // create filename: "rst/rank_YYYYYYY/file_basename" + "." + XXXXX + ".rst"
    // where YYYYYYY = 8-digit rank number
    // where XXXXX = 5-digit file_number
    char rank_dir[20];
    char number[7];
    std::snprintf(number, sizeof(number), ".%05d", out_params.file_number);
    std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d/", global_variable::my_rank);
    fname = std::string("rst/") + std::string(rank_dir) + out_params.file_basename
      + number + ".rst";

    // Debugging output to check directory and filename
    // std::cout << "Rank " << global_variable::my_rank << " generated filename: "
    //           << fname << std::endl;
  } else {
    // Existing behavior: single restart file
    // create filename: "rst/file_basename" + "." + XXXXX + ".rst"
    // where XXXXX = 5-digit file_number
    char number[7];
    std::snprintf(number, sizeof(number), ".%05d", out_params.file_number);
    fname = std::string("rst/") + out_params.file_basename + number + ".rst";
  }
  // increment counters now so values for *next* dump are stored in restart file
  out_params.file_number++;
  // Driver::Finalize's restart file before tlim keeps the schedule
  // (the next link writes its next restart file when the continuous run would)
  if (global_variable::final_outputs && out_params.last_time >= 0.0) {
    // last_time unchanged
  } else if (out_params.last_time < 0.0) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  pin->SetInteger(out_params.block_name, "file_number", out_params.file_number);
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);

  // create string holding input parameters (copy of input file)
  std::stringstream ost;
  pin->ParameterDump(ost);
  std::string sbuf = ost.str();

  //--- STEP 1.  Root process writes header data (input file, critical variables)
  // Input file data is read by ParameterInput on restart, and the remaining header
  // variables are read in Mesh::BuildTreeFromRestart()

  // open file and  write the header; this part is serial
  IOWrapper resfile;
  resfile.Open(fname.c_str(), IOWrapper::FileMode::write, single_file_per_rank);
  // rank 0's progress through the dump, and its watchdog: every <job>/restart_warn_s
  // seconds (default 600; 0 off) without a progress line it warns, naming the last one,
  // so a stalled file system is never silent
  const auto t0w = std::chrono::steady_clock::now();
  auto sec = [&t0w]() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0w).count();
  };
  const double warn_s = pin->GetOrAddReal("job", "restart_warn_s", 600.0);
  std::mutex wmx;
  std::condition_variable wcv;
  bool wdone = false;
  double tlast = 0.0;
  std::string wlast = "the open";
  std::thread wd;
  if (global_variable::my_rank == 0 && warn_s > 0.0) {
    wd = std::thread([&]() {
      std::unique_lock<std::mutex> lk(wmx);
      while (!wdone) {
        if (wcv.wait_for(lk, std::chrono::duration<double>(warn_s)) ==
            std::cv_status::timeout && !wdone && sec() - tlast >= warn_s) {
          std::cout << "CLEAN WARNING restart dump " << fname << ": no progress for "
                    << static_cast<int>(sec() - tlast) << " s after " << wlast
                    << std::endl;
        }
      }
    });
  }
  auto prog = [&](const char *what) {
    if (global_variable::my_rank != 0) return;
    std::lock_guard<std::mutex> lk(wmx);
    tlast = sec();
    wlast = what;
    std::cout << "<restart> " << fname << ": " << what << " written at " << std::fixed
              << std::setprecision(1) << tlast << " s" << std::endl;
    std::cout.unsetf(std::ios_base::floatfield);
  };
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    // output the input parameters (input file)
    resfile.Write_any_type(sbuf.c_str(), sbuf.size(), "byte", single_file_per_rank);

    // output Mesh information
    resfile.Write_any_type(&(pm->nmb_total), (sizeof(int)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->root_level), (sizeof(int)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->mesh_size), (sizeof(RegionSize)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->mesh_indcs), (sizeof(RegionIndcs)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->mb_indcs), (sizeof(RegionIndcs)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->time), (sizeof(Real)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->dt), (sizeof(Real)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->ncycle), (sizeof(int)), "byte",
                            single_file_per_rank);
  }
  //--- STEP 2.  Root process writes list of logical locations and cost of MeshBlocks
  // This data read in Mesh::BuildTreeFromRestart()

  if (global_variable::my_rank == 0 || single_file_per_rank) {
    resfile.Write_any_type(&(pm->lloc_eachmb[0]),(pm->nmb_total)*sizeof(LogicalLocation),
                           "byte", single_file_per_rank);
    resfile.Write_any_type(&(pm->cost_eachmb[0]), (pm->nmb_total)*sizeof(float),
                           "byte", single_file_per_rank);
  }

  //--- STEP 3.  Root process writes internal state of objects that require it
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    // store z4c information
    if (pz4c != nullptr) {
      resfile.Write_any_type(&(pz4c->last_output_time), sizeof(Real), "byte",
                             single_file_per_rank);
    }
    // output puncture tracker data
    if (nco > 0) {
      for (auto & pt : pz4c->ptracker) {
        resfile.Write_any_type(pt->GetPos(), 3*sizeof(Real), "byte",
                               single_file_per_rank);
      }
    }
    // turbulence driver internal RNG
    if (pturb != nullptr) {
      resfile.Write_any_type(&(pturb->rstate), sizeof(RNG_State), "byte",
                             single_file_per_rank);
    }
    // the pgen's own relaxed state, marker then length then payload
    if (!pgen_state.empty()) {
      IOWrapperSizeT nb = pgen_state.size();
      resfile.Write_any_type(&(kPgenRstMagic[0]), sizeof(kPgenRstMagic), "byte",
                             single_file_per_rank);
      resfile.Write_any_type(&nb, sizeof(IOWrapperSizeT), "byte", single_file_per_rank);
      resfile.Write_any_type(pgen_state.data(), nb, "byte", single_file_per_rank);
    }
    // the mode-3 warm-start header, marker then length then payload, read by the same
    // eight-byte peek the pgen state block uses (utils/two_stream_warm_rst.hpp)
    if (nwarm > 0) {
      IOWrapperSizeT nb = sizeof(warm_hdr);
      resfile.Write_any_type(&(two_stream_rt::kRtWarmRstMagic[0]),
                             sizeof(two_stream_rt::kRtWarmRstMagic), "byte",
                             single_file_per_rank);
      resfile.Write_any_type(&nb, sizeof(IOWrapperSizeT), "byte", single_file_per_rank);
      resfile.Write_any_type(&(warm_hdr[0]), nb, "byte", single_file_per_rank);
    }
    // the implicit_predictor header, same marked form, behind the warm-start one
    if (npred > 0) {
      IOWrapperSizeT nb = sizeof(pred_hdr);
      resfile.Write_any_type(&(radvet::kVetPredRstMagic[0]),
                             sizeof(radvet::kVetPredRstMagic), "byte",
                             single_file_per_rank);
      resfile.Write_any_type(&nb, sizeof(IOWrapperSizeT), "byte", single_file_per_rank);
      resfile.Write_any_type(&(pred_hdr[0]), nb, "byte", single_file_per_rank);
    }
    // the hesdirk2 slope header, same marked form, behind the predictor one
    if (nt2 > 0) {
      IOWrapperSizeT nb = sizeof(t2_hdr);
      resfile.Write_any_type(&(radvet::kVetTime2RstMagic[0]),
                             sizeof(radvet::kVetTime2RstMagic), "byte",
                             single_file_per_rank);
      resfile.Write_any_type(&nb, sizeof(IOWrapperSizeT), "byte", single_file_per_rank);
      resfile.Write_any_type(&(t2_hdr[0]), nb, "byte", single_file_per_rank);
    }
    if (wonep) {
      IOWrapperSizeT nb = onep_nb;
      resfile.Write_any_type(&(radvet::kVetOnePassRstMagic[0]),
                             sizeof(radvet::kVetOnePassRstMagic), "byte",
                             single_file_per_rank);
      resfile.Write_any_type(&nb, sizeof(IOWrapperSizeT), "byte", single_file_per_rank);
      resfile.Write_any_type(&(onep_hdr[0]), nb, "byte", single_file_per_rank);
    }
    if (wabe) {
      IOWrapperSizeT nb = sizeof(abe_hdr);
      resfile.Write_any_type(&(radvet::kVetAutoBeRstMagic[0]),
                             sizeof(radvet::kVetAutoBeRstMagic), "byte",
                             single_file_per_rank);
      resfile.Write_any_type(&nb, sizeof(IOWrapperSizeT), "byte", single_file_per_rank);
      resfile.Write_any_type(&(abe_hdr[0]), nb, "byte", single_file_per_rank);
    }
    // the internal energy header, same marked form, behind all the others
    if (neint > 0) {
      IOWrapperSizeT nb = sizeof(eint_hdr);
      resfile.Write_any_type(&(kEintRstMagic[0]), sizeof(kEintRstMagic), "byte",
                             single_file_per_rank);
      resfile.Write_any_type(&nb, sizeof(IOWrapperSizeT), "byte", single_file_per_rank);
      resfile.Write_any_type(&(eint_hdr[0]), nb, "byte", single_file_per_rank);
    }
    // the implicit-ck state header, same marked form, behind all the others
    if (nck > 0) {
      IOWrapperSizeT nb = ck_rst_hdr_bytes.size();
      resfile.Write_any_type(&(two_stream_rt::kCkRstMagic[0]),
                             sizeof(two_stream_rt::kCkRstMagic), "byte",
                             single_file_per_rank);
      resfile.Write_any_type(&nb, sizeof(IOWrapperSizeT), "byte", single_file_per_rank);
      resfile.Write_any_type(ck_rst_hdr_bytes.data(), nb, "byte", single_file_per_rank);
    }
    // the thin_relax closure header, same marked form, behind the implicit-ck one
    if (nctr > 0) {
      IOWrapperSizeT nb = sizeof(ctr_hdr);
      resfile.Write_any_type(&(radvet::kVetCtrRstMagic[0]),
                             sizeof(radvet::kVetCtrRstMagic),
                             "byte", single_file_per_rank);
      resfile.Write_any_type(&nb, sizeof(IOWrapperSizeT), "byte", single_file_per_rank);
      resfile.Write_any_type(&(ctr_hdr[0]), nb, "byte", single_file_per_rank);
    }
    // the signal-speed header, same marked form, behind the thin_relax one
    if (nrss > 0) {
      IOWrapperSizeT nb = sizeof(rss_hdr);
      resfile.Write_any_type(&(radvet::kVetRssRstMagic[0]),
                             sizeof(radvet::kVetRssRstMagic),
                             "byte", single_file_per_rank);
      resfile.Write_any_type(&nb, sizeof(IOWrapperSizeT), "byte", single_file_per_rank);
      resfile.Write_any_type(&(rss_hdr[0]), nb, "byte", single_file_per_rank);
    }
    // the conservative cell-F slope header, same marked form, behind the signal-speed one
    if (wcsk) {
      IOWrapperSizeT nb = sizeof(csk_hdr);
      resfile.Write_any_type(&(radvet::kVetCskRstMagic[0]),
                             sizeof(radvet::kVetCskRstMagic),
                             "byte", single_file_per_rank);
      resfile.Write_any_type(&nb, sizeof(IOWrapperSizeT), "byte", single_file_per_rank);
      resfile.Write_any_type(&(csk_hdr[0]), nb, "byte", single_file_per_rank);
    }
    // the vet_rotate_avg header, same marked form, behind the csk one
    if (nra > 0) {
      IOWrapperSizeT nb = sizeof(ra_hdr);
      resfile.Write_any_type(&(radvet::kVetRaRstMagic[0]), sizeof(radvet::kVetRaRstMagic),
                             "byte", single_file_per_rank);
      resfile.Write_any_type(&nb, sizeof(IOWrapperSizeT), "byte", single_file_per_rank);
      resfile.Write_any_type(&(ra_hdr[0]), nb, "byte", single_file_per_rank);
    }
    // the am_vector spin header, same marked form, behind the vet_rotate_avg one
    if (nams > 0) {
      IOWrapperSizeT nb = sizeof(ams_hdr);
      resfile.Write_any_type(&(hydro::kAmSpinRstMagic[0]), sizeof(hydro::kAmSpinRstMagic),
                             "byte", single_file_per_rank);
      resfile.Write_any_type(&nb, sizeof(IOWrapperSizeT), "byte", single_file_per_rank);
      resfile.Write_any_type(&(ams_hdr[0]), nb, "byte", single_file_per_rank);
    }
    // the per-cell pgen slab header, same marked form, behind the spin one
    if (npgs > 0) {
      IOWrapperSizeT nb = sizeof(pgs_hdr);
      resfile.Write_any_type(&(kPgenSlabRstMagic[0]), sizeof(kPgenSlabRstMagic),
                             "byte", single_file_per_rank);
      resfile.Write_any_type(&nb, sizeof(IOWrapperSizeT), "byte", single_file_per_rank);
      resfile.Write_any_type(&(pgs_hdr[0]), nb, "byte", single_file_per_rank);
    }
    // the rad_vet applied face flux header, same marked form, behind the pgen slab one
    if (nfax > 0) {
      IOWrapperSizeT nb = sizeof(fax_hdr);
      resfile.Write_any_type(&(radvet::kVetFaxRstMagic[0]),
                             sizeof(radvet::kVetFaxRstMagic),
                             "byte", single_file_per_rank);
      resfile.Write_any_type(&nb, sizeof(IOWrapperSizeT), "byte", single_file_per_rank);
      resfile.Write_any_type(&(fax_hdr[0]), nb, "byte", single_file_per_rank);
    }
    // the vet_x1_periodic warm-start header, behind the fax one
    if (nx1e > 0) {
      IOWrapperSizeT nb = sizeof(x1e_hdr);
      resfile.Write_any_type(&(radvet::kVetX1eRstMagic[0]),
                             sizeof(radvet::kVetX1eRstMagic),
                             "byte", single_file_per_rank);
      resfile.Write_any_type(&nb, sizeof(IOWrapperSizeT), "byte", single_file_per_rank);
      resfile.Write_any_type(&(x1e_hdr[0]), nb, "byte", single_file_per_rank);
    }
  }

  //--- STEP 4.  All ranks write data over all MeshBlocks (5D arrays) in parallel
  // This data read in ProblemGenerator constructor for restarts

  // total size of all cell-centered variables and face-centered fields to be written by
  // this rank
  // the general-EOS temperature cache is appended at the END of each MeshBlock record,
  // so that a file written without it (any run with an ideal gas, and every file written
  // before this was added) differs only in its tail and is still readable -- see the
  // size check in pgen.cpp's restart constructor.
  bool wt_hyd = (phydro != nullptr) && phydro->peos->eos_data.IsGeneral();
  bool wt_mhd = (pmhd != nullptr) && pmhd->peos->eos_data.IsGeneral();
  IOWrapperSizeT data_size = 0;
  // bytes of one cell-centred variable over one MeshBlock
  const IOWrapperSizeT cell_bytes =
      static_cast<IOWrapperSizeT>(nout1)*nout2*nout3*sizeof(Real);
  if (phydro != nullptr) {
    data_size += nout1*nout2*nout3*nhydro*sizeof(Real); // hydro u0
  }
  if (pmhd != nullptr) {
    data_size += nout1*nout2*nout3*nmhd*sizeof(Real);   // mhd u0
    data_size += (nout1+1)*nout2*nout3*sizeof(Real);    // mhd b0.x1f
    data_size += nout1*(nout2+1)*nout3*sizeof(Real);    // mhd b0.x2f
    data_size += nout1*nout2*(nout3+1)*sizeof(Real);    // mhd b0.x3f
  }
  if (prad != nullptr) {
    data_size += nout1*nout2*nout3*nrad*sizeof(Real);   // radiation i0
  }
  if (pradvet != nullptr) {
    data_size += nout1*nout2*nout3*nm1*sizeof(Real);    // rad_vet u0
    {
      data_size += (nout1+1)*nout2*nout3*sizeof(Real);  // rad_vet f0x1 (milestone 3a)
    }
    if (pradvet->trans_on) {
      data_size += nout1*(nout2+1)*nout3*sizeof(Real);  // rad_vet f0x2 (3b phase B)
      if (pradvet->trans_x3) {
        data_size += nout1*nout2*(nout3+1)*sizeof(Real);  // rad_vet f0x3
      }
    }
  }
  if (pturb != nullptr) {
    data_size += nout1*nout2*nout3*nforce*sizeof(Real); // forcing
  }
  if (pz4c != nullptr) {
    data_size += nout1*nout2*nout3*nz4c*sizeof(Real);   // z4c u0
  } else if (padm != nullptr) {
    data_size += nout1*nout2*nout3*nadm*sizeof(Real);   // adm u_adm
  }
  if (wt_hyd) {
    data_size += nout1*nout2*nout3*sizeof(Real);        // hydro wtemp
  }
  if (wt_mhd) {
    data_size += nout1*nout2*nout3*sizeof(Real);        // mhd wtemp
  }
  if (wt_hyd) {
    data_size += 2*nout1*nout2*nout3*sizeof(Real);      // hydro wder (IDPR, IDG1)
  }
  if (wt_mhd) {
    data_size += 2*nout1*nout2*nout3*sizeof(Real);      // mhd wder (IDPR, IDG1)
  }
  if (nwarm > 0) {
    data_size += nwarm*cell_bytes;   // rt_c3bp (+ rt_c3bp2)
  }
  if (npred > 0) {
    data_size += npred*cell_bytes;   // rad_vet ipred
  }
  if (nt2 > 0) {
    data_size += nt2*cell_bytes;   // rad_vet hesdirk2 slope
  }
  // hydro / mhd w0(IEN), + the hydro w0 density and velocities
  data_size += (neint + (we_hyd ? 4 : 0))*cell_bytes;
  data_size += nck*cell_bytes;   // implicit-ck state
  data_size += nctr*cell_bytes;   // rad_vet thin_relax ctr_mem
  data_size += nrss*cell_bytes;   // rad_signal_speed inputs
  data_size += ncsk*cell_bytes;   // conservative csk1
  data_size += nra*cell_bytes;   // vet_rotate_avg history
  data_size += nams*cell_bytes;   // am_vector spin (+ amk1)
  data_size += npgs*cell_bytes;   // per-cell pgen slabs
  data_size += nfax*cell_bytes;   // fax*
  data_size += nx1e*cell_bytes;   // vet_x1_periodic warm start
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    resfile.Write_any_type(&(data_size), sizeof(IOWrapperSizeT), "byte",
                            single_file_per_rank);
  }

  // calculate size of data written in Steps 1-2 above
  IOWrapperSizeT step1size = sbuf.size()*sizeof(char) + 3*sizeof(int) + 2*sizeof(Real) +
                             sizeof(RegionSize) + 2*sizeof(RegionIndcs);
  IOWrapperSizeT step2size = (pm->nmb_total)*(sizeof(LogicalLocation) + sizeof(float));

  IOWrapperSizeT step3size = 3*nco*sizeof(Real);
  if (pz4c != nullptr) step3size += sizeof(Real);
  if (pturb != nullptr) step3size += sizeof(RNG_State);
  if (!pgen_state.empty()) {
    step3size += sizeof(kPgenRstMagic) + sizeof(IOWrapperSizeT) + pgen_state.size();
  }
  if (nwarm > 0) {
    step3size += sizeof(two_stream_rt::kRtWarmRstMagic) + sizeof(IOWrapperSizeT)
                 + sizeof(warm_hdr);
  }
  if (npred > 0) {
    step3size += sizeof(radvet::kVetPredRstMagic) + sizeof(IOWrapperSizeT)
                 + sizeof(pred_hdr);
  }
  if (nt2 > 0) {
    step3size += sizeof(radvet::kVetTime2RstMagic) + sizeof(IOWrapperSizeT)
                 + sizeof(t2_hdr);
  }
  if (wonep) {
    step3size += sizeof(radvet::kVetOnePassRstMagic) + sizeof(IOWrapperSizeT)
                 + onep_nb;
  }
  if (wabe) {
    step3size += sizeof(radvet::kVetAutoBeRstMagic) + sizeof(IOWrapperSizeT)
                 + sizeof(abe_hdr);
  }
  if (neint > 0) {
    step3size += sizeof(kEintRstMagic) + sizeof(IOWrapperSizeT) + sizeof(eint_hdr);
  }
  if (nck > 0) {
    step3size += sizeof(two_stream_rt::kCkRstMagic) + sizeof(IOWrapperSizeT)
                 + ck_rst_hdr_bytes.size();
  }
  if (nctr > 0) {
    step3size += sizeof(radvet::kVetCtrRstMagic) + sizeof(IOWrapperSizeT) +
        sizeof(ctr_hdr);
  }
  if (nrss > 0) {
    step3size += sizeof(radvet::kVetRssRstMagic) + sizeof(IOWrapperSizeT) +
        sizeof(rss_hdr);
  }
  if (wcsk) {
    step3size += sizeof(radvet::kVetCskRstMagic) + sizeof(IOWrapperSizeT) +
        sizeof(csk_hdr);
  }
  if (nra > 0) {
    step3size += sizeof(radvet::kVetRaRstMagic) + sizeof(IOWrapperSizeT) + sizeof(ra_hdr);
  }
  if (nams > 0) {
    step3size += sizeof(hydro::kAmSpinRstMagic) + sizeof(IOWrapperSizeT) +
                 sizeof(ams_hdr);
  }
  if (npgs > 0) {
    step3size += sizeof(kPgenSlabRstMagic) + sizeof(IOWrapperSizeT) + sizeof(pgs_hdr);
  }
  if (nfax > 0) {
    step3size += sizeof(radvet::kVetFaxRstMagic) + sizeof(IOWrapperSizeT) +
        sizeof(fax_hdr);
  }
  if (nx1e > 0) {
    step3size += sizeof(radvet::kVetX1eRstMagic) + sizeof(IOWrapperSizeT) +
        sizeof(x1e_hdr);
  }

  // write cell-centered variables in parallel
  IOWrapperSizeT offset_myrank = (step1size + step2size + step3size
                                  + sizeof(IOWrapperSizeT));

  if (!single_file_per_rank) {
    offset_myrank += data_size*(pm->gids_eachrank[global_variable::my_rank]);
  }

  IOWrapperSizeT myoffset = offset_myrank;

  // write cell-centered variables, one MeshBlock at a time (but parallelized over all
  // ranks). MeshBlocks are written seperately to reduce number of data elements per write
  // call, to avoid exceeding 2^31 limit for very large grids per MPI rank.
  if (phydro != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_hyd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered hydro data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_hyd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered hydro data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nhydro*sizeof(Real); // hydro u0
    myoffset = offset_myrank;
  }
  if (pmhd != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_mhd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered mhd data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_mhd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered mhd data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nmhd*sizeof(Real);   // mhd u0
    myoffset = offset_myrank;

    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to x1-face field
        auto x1fptr = Kokkos::subview(outfield.x1f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        int fldcnt = x1fptr.size();
        if (resfile.Write_any_type_at_all(x1fptr.data(),fldcnt,myoffset,"Real",
                                          single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x1f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        // get ptr to x2-face field
        auto x2fptr = Kokkos::subview(outfield.x2f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = x2fptr.size();
        if (resfile.Write_any_type_at_all(x2fptr.data(),fldcnt,myoffset,"Real",
                                          single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x2f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        // get ptr to x3-face field
        auto x3fptr = Kokkos::subview(outfield.x3f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = x3fptr.size();
        if (resfile.Write_any_type_at_all(x3fptr.data(),fldcnt,myoffset,"Real",
                                          single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x3f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        myoffset += data_size-(x1fptr.size()+x2fptr.size()+x3fptr.size())*sizeof(Real);

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to x1-face field
        auto x1fptr = Kokkos::subview(outfield.x1f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        int fldcnt = x1fptr.size();
        if (resfile.Write_any_type_at(x1fptr.data(),fldcnt,myoffset,"Real",
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x1f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        // get ptr to x2-face field
        auto x2fptr = Kokkos::subview(outfield.x2f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = x2fptr.size();
        if (resfile.Write_any_type_at(x2fptr.data(),fldcnt,myoffset,"Real",
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x2f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        // get ptr to x3-face field
        auto x3fptr = Kokkos::subview(outfield.x3f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = x3fptr.size();
        if (resfile.Write_any_type_at(x3fptr.data(),fldcnt,myoffset,"Real",
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x3f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        myoffset += data_size-(x1fptr.size()+x2fptr.size()+x3fptr.size())*sizeof(Real);
      }
    }
    offset_myrank += (nout1+1)*nout2*nout3*sizeof(Real);    // mhd b0.x1f
    offset_myrank += nout1*(nout2+1)*nout3*sizeof(Real);    // mhd b0.x2f
    offset_myrank += nout1*nout2*(nout3+1)*sizeof(Real);    // mhd b0.x3f
    myoffset = offset_myrank;
  }

  if (prad != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_rad, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered rad data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_rad, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(),mbcnt,myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered rad data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nrad*sizeof(Real);   // radiation i0
    myoffset = offset_myrank;
  }

  if (pradvet != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        auto mbptr = Kokkos::subview(outarray_m1, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered rad_vet data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        auto mbptr = Kokkos::subview(outarray_m1, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(),mbcnt,myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered rad_vet data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nm1*sizeof(Real);    // rad_vet u0
    myoffset = offset_myrank;

    // milestone 3a / 3b phase B: the persistent face fluxes, same loop, same order.
    // ONE writer for the three face arrays, so that the order here cannot drift from the
    // order of the reader in src/pgen/pgen.cpp.
    auto wrface = [&](HostArray4D<Real> &arr) {
      for (int m=0;  m<noutmbs_max; ++m) {
        if (m < noutmbs_min) {
          auto mbptr = Kokkos::subview(arr, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
          int mbcnt = mbptr.size();
          if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                            single_file_per_rank) != mbcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "rad_vet face data not written correctly to rst file, "
            << "restart file is broken." << std::endl;
            exit(EXIT_FAILURE);
          }
          myoffset += data_size;
        } else if (m < pm->nmb_thisrank) {
          auto mbptr = Kokkos::subview(arr, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
          int mbcnt = mbptr.size();
          if (resfile.Write_any_type_at(mbptr.data(),mbcnt,myoffset,"Real",
                                        single_file_per_rank) != mbcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "rad_vet face data not written correctly"
                      << " to rst file, restart file is broken." << std::endl;
            exit(EXIT_FAILURE);
          }
          myoffset += data_size;
        }
      }
    };
    {
      for (int m=0;  m<noutmbs_max; ++m) {
        if (m < noutmbs_min) {
          auto mbptr = Kokkos::subview(outarray_m1f, m, Kokkos::ALL, Kokkos::ALL,
                                       Kokkos::ALL);
          int mbcnt = mbptr.size();
          if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                            single_file_per_rank) != mbcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "rad_vet face data not written correctly to rst file, "
            << "restart file is broken." << std::endl;
            exit(EXIT_FAILURE);
          }
          myoffset += data_size;
        } else if (m < pm->nmb_thisrank) {
          auto mbptr = Kokkos::subview(outarray_m1f, m, Kokkos::ALL, Kokkos::ALL,
                                       Kokkos::ALL);
          int mbcnt = mbptr.size();
          if (resfile.Write_any_type_at(mbptr.data(),mbcnt,myoffset,"Real",
                                        single_file_per_rank) != mbcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "rad_vet face data not written correctly"
                      << " to rst file, restart file is broken." << std::endl;
            exit(EXIT_FAILURE);
          }
          myoffset += data_size;
        }
      }
      offset_myrank += (nout1+1)*nout2*nout3*sizeof(Real);   // rad_vet f0x1
      myoffset = offset_myrank;
    }
    if (pradvet->trans_on) {
      wrface(outarray_m1f2);
      offset_myrank += nout1*(nout2+1)*nout3*sizeof(Real);   // rad_vet f0x2
      myoffset = offset_myrank;
      if (pradvet->trans_x3) {
        wrface(outarray_m1f3);
        offset_myrank += nout1*nout2*(nout3+1)*sizeof(Real);   // rad_vet f0x3
        myoffset = offset_myrank;
      }
    }
  }

  if (pturb != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_force, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered turb data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_force, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered turb data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nforce*sizeof(Real); // forcing
    myoffset = offset_myrank;
  }

  if (pz4c != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_z4c, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered z4c data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_z4c, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered z4c data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nz4c*sizeof(Real); // z4c u0
    myoffset = offset_myrank;
  } else if (padm != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_adm, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered adm data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_adm, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered adm data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nadm*sizeof(Real); // adm u_adm
    myoffset = offset_myrank;
  }

  // write the general-EOS temperature cache last, one MeshBlock at a time, exactly as
  // the arrays above.  Same loop for hydro and MHD, so it is written once here.
  auto write_wtemp = [&](const auto &a, const char *what) {
    for (int m=0;  m<noutmbs_max; ++m) {
      if (m < noutmbs_min) {
        auto mbptr = Kokkos::subview(a, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << what << " cache not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      } else if (m < pm->nmb_thisrank) {
        auto mbptr = Kokkos::subview(a, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << what << " cache not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*sizeof(Real);
    myoffset = offset_myrank;
  };
  // a group of n slabs (a(m, n, k, j, i)): with independent I/O a MeshBlock's slabs are
  // contiguous in memory and in its record, so they go out in one write per MeshBlock (in
  // pieces of at most 1 GiB: ROMIO fails a single request of 2^31 bytes); same bytes as the
  // slab by slab loop
  using PackFn = std::function<void(DvceArray5D<Real> &, int, int)>;
  auto write_group = [&](const auto &a, const int n, const char *what,
                         const PackFn &pk = PackFn()) {
    if (n <= 0) return;
    if (IOWrapper::collective) {
      for (int c=0; c<n; ++c) {
        write_wtemp(Kokkos::subview(a, Kokkos::ALL, c, Kokkos::ALL, Kokkos::ALL,
                                    Kokkos::ALL), what);
      }
      prog(what);
      return;
    }
    const IOWrapperSizeT cnt = static_cast<IOWrapperSizeT>(n)*nout1*nout2*nout3;
    // (streamed, the time2 record: packed MeshBlock by MeshBlock in windows of nw slabs)
    const bool st = (a.extent_int(0) < pm->nmb_thisrank);
    const int nw = st ? a.extent_int(1) : n;
    const IOWrapperSizeT ns = static_cast<IOWrapperSizeT>(nout1)*nout2*nout3;
    for (int m=0; m<pm->nmb_thisrank; ++m) {
      for (int w0=0; w0<n; w0+=nw) {
        if (st) {
          DvceArray5D<Real> t(const_cast<Real *>(a.data()), 1, nw, nout3, nout2, nout1);
          pk(t, m, w0);
          Kokkos::fence();
        }
        const Real *p = Kokkos::subview(a, st ? 0 : m, Kokkos::ALL, Kokkos::ALL,
                                        Kokkos::ALL, Kokkos::ALL).data();
        const IOWrapperSizeT cw = std::min(nw, n - w0)*ns, o = w0*ns;
        for (IOWrapperSizeT q=0; q<cw; q+=(IOWrapperSizeT{1} << 27)) {
          const IOWrapperSizeT c = std::min(cw - q, IOWrapperSizeT{1} << 27);
          if (resfile.Write_any_type_at(p + q, c, offset_myrank + m*data_size +
                                        (o + q)*sizeof(Real), "Real",
                                        single_file_per_rank) != c) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << what << " not written correctly to rst file, "
                      << "restart file is broken." << std::endl;
            exit(EXIT_FAILURE);
          }
        }
      }
    }
    offset_myrank += cnt*sizeof(Real);
    myoffset = offset_myrank;
    prog(what);
  };
  prog("the cell and face variables");
  if (wt_hyd) { write_wtemp(outarray_wth, "hydro"); }
  if (wt_mhd) { write_wtemp(outarray_wtm, "mhd"); }
  // the derived cache follows it, same layout, same loop (see outarray_wdp)
  if (wt_hyd) {
    write_wtemp(outarray_wdp, "hydro pressure");
    write_wtemp(outarray_wdg, "hydro Gamma_1");
  }
  if (wt_mhd) {
    write_wtemp(outarray_wdpm, "mhd pressure");
    write_wtemp(outarray_wdgm, "mhd Gamma_1");
  }
  // and the warm-start history behind that, same layout, same loop
  if (nwarm > 0) { write_wtemp(outarray_wm1, "rt warm start"); }
  if (nwarm > 1) { write_wtemp(outarray_wm2, "rt warm start 2"); }
  // and the implicit_predictor increment behind that, one slab per channel
  write_group(outarray_m1p, npred, "rad_vet predictor");
  // and the hesdirk2 slope behind that
  write_group(outarray_m1t, nt2, "rad_vet hesdirk2 slope",
              [&](DvceArray5D<Real> &t, int m, int w0) {
                pradvet->Time2RstPack(t, 1, m, w0);});
  // and the internal energy last of all (kEintRstMagic)
  if (we_hyd) { write_wtemp(outarray_weh, "hydro internal energy"); }
  // (w0's IDN, IVX, IVY, IVZ are its first 4 channels: contiguous per MeshBlock)
  write_group(outarray_wpv, we_hyd ? 4 : 0, "hydro w0 density and velocities");
  if (we_mhd) { write_wtemp(outarray_wem, "mhd internal energy"); }
  // and the implicit-ck state behind it (utils/two_stream_ck_rst.hpp)
  for (int n=0; n<nck; ++n) {
    write_wtemp(Kokkos::subview(outarray_ck, n, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                Kokkos::ALL), "ck restart state");
  }
  // and the thin_relax closure memory last of all (kVetCtrRstMagic)
  write_group(outarray_m1c, nctr, "rad_vet thin_relax closure");
  // and the rad_signal_speed inputs behind it (kVetRssRstMagic)
  write_group(outarray_m1s, nrss, "rad_signal_speed inputs",
              [&](DvceArray5D<Real> &t, int m, int) {
                pradvet->RssRstPack(t, 1, pradvet->RssRstMode(), m);});
  // and the vet_gd_async pending source behind it (kVetVgdRstMagic)
  // and the conservative cell-F slope behind it (kVetCskRstMagic)
  write_group(outarray_m1k, ncsk, "rad_vet conservative cell-F slope");
  // and the vet_rotate_avg history behind it (kVetRaRstMagic)
  write_group(outarray_m1r, nra, "vet_rotate_avg history");
  // and the am_vector spin (+ amk1) behind it (hydro::kAmSpinRstMagic)
  write_group(outarray_ams, nams, "am_vector spin");
  // and the per-cell pgen slabs behind it (kPgenSlabRstMagic)
  write_group(outarray_pgs, npgs, "pgen slab");
  // and the rad_vet applied face fluxes last of all (kVetFaxRstMagic)
  write_group(outarray_fax, nfax, "rad_vet applied face flux");
  // and the vet_x1_periodic warm-start plane behind them (kVetX1eRstMagic)
  write_group(outarray_x1e, nx1e, "rad_vet vet_x1 warm start");

  // close file, clean up
  resfile.Close(single_file_per_rank);
  prog("all (file closed)");
  if (wd.joinable()) {
    { std::lock_guard<std::mutex> lk(wmx); wdone = true; }
    wcv.notify_all();
    wd.join();
  }
  // drop the host staging copies (filled again by the next LoadOutputData,
  // whose Kokkos::realloc allocates anew anyway): on a unified-memory APU (MI300A) they
  // held a full second copy of the restart data in HBM for the rest of the run
  outarray = HostArray5D<Real>();
  outarray_hyd = outarray_mhd = outarray_rad = outarray_m1 = outarray_force =
      outarray_z4c = outarray_adm = outarray_ck = outarray_m1p = outarray_m1t =
      outarray_m1c = outarray_m1s = outarray_m1r = outarray_m1k =
      outarray_ams = outarray_pgs = outarray_fax = outarray_x1e = outarray_wpv =
      HostArray5D<Real>();
  outarray_m1t_d = DvceArray5D<Real>();
  rst_keep_.clear();
  outarray_m1f = outarray_m1f2 = outarray_m1f3 = outarray_wth = outarray_wtm =
      outarray_wdp = outarray_wdg = outarray_wdpm = outarray_wdgm = outarray_weh =
      outarray_wem = outarray_wm1 = outarray_wm2 = HostArray4D<Real>();
  outfield.x1f = outfield.x2f = outfield.x3f = HostArray4D<Real>();


  return;
}
