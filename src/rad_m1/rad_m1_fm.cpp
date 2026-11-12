//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_fm.cpp
//! \brief fmflux-1113 <rad_m1>/implicit_hr_cmem_fm: the state of the switch-free
//! momentum face flux F_m (rad_m1.hpp).  F_m is evolved by the face kernels of the
//! implicit solve (rad_m1_implicit.cpp: x1 m1_impl_face, x2 / x3 in
//! ImplicitTransverseTerms) as F_m' = F_a' + coef (F_m^n - F_c^n); this file keeps the
//! rest of its life cycle in step with the F_c memory f0x*: the solve-start copy (+ the
//! cubed-sphere seam average), the hesdirk2 stage copies, Heun average, increments and
//! slopes, the posfix drop and the rotation correction's K1.  Nothing here runs unless
//! impl_cm_fm is set.

#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

namespace radm1 {

namespace {
void FmAlloc4(DvceArray4D<Real> &a, const DvceArray4D<Real> &ref) {
  if (ref.extent_int(0) == 0) return;
  if (a.extent(0) != ref.extent(0) || a.extent(1) != ref.extent(1) ||
      a.extent(2) != ref.extent(2) || a.extent(3) != ref.extent(3)) {
    Kokkos::realloc(a, ref.extent(0), ref.extent(1), ref.extent(2), ref.extent(3));
    Kokkos::deep_copy(a, 0.0);
  }
}
void FmAlloc5(DvceArray5D<Real> &a, const DvceArray5D<Real> &ref) {
  if (ref.extent_int(0) == 0) return;
  if (a.extent(0) != ref.extent(0) || a.extent_int(1) != 3 ||
      a.extent(2) != ref.extent(2) || a.extent(3) != ref.extent(3) ||
      a.extent(4) != ref.extent(4)) {
    Kokkos::realloc(a, ref.extent(0), 3, ref.extent(2), ref.extent(3), ref.extent(4));
    Kokkos::deep_copy(a, 0.0);
  }
}
void FmCopy(DvceArray4D<Real> &dst, const DvceArray4D<Real> &src) {
  if (dst.extent_int(0) > 0 && dst.extent(0) == src.extent(0) &&
      dst.extent(3) == src.extent(3)) {
    Kokkos::deep_copy(DevExeSpace(), dst, src);
  }
}
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::FmEnsure
//! \brief allocate F_m's arrays (as f0x*, t2f*, t2k1) and, once, start F_m from the
//! applied fluxes fax* (from f0x* while implicit_hr_cmem_fainit has not run yet: the
//! value it will give fax*) and its K1 from the F channels of t2k1 (a restart file
//! carries neither: F_m^n - F_c^n continues as F_a^n - F_c^n)

void RadiationM1::FmEnsure() {
  if (!impl_cm_fm) return;
  if (mr_every > 1 || t2_k1reset) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<rad_m1>/implicit_hr_cmem_fm: not with implicit_mr_every > 1 or "
              << "time2_k1_reset (their face-slope edits are not carried to F_m)"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  FmAlloc4(fm1, f0x1);
  FmAlloc4(fm1n, f0x1);
  FmAlloc4(t2fm1, t2f1);
  if (trans_on) {
    FmAlloc4(fm2, f0x2);
    FmAlloc4(fm2n, f0x2);
    FmAlloc4(t2fm2, t2f2);
    if (trans_x3) {
      FmAlloc4(fm3, f0x3);
      FmAlloc4(fm3n, f0x3);
      FmAlloc4(t2fm3, t2f3);
    }
  }
  FmAlloc5(t2km1, t2k1);
  FmAlloc5(t2km2, t2k1);
  FmAlloc5(t2incm, t2k1);
  if (fm_init_done) return;
  fm_init_done = true;
  const bool fromf = impl_cm_fainit && !impl_cm_fainit_done;
  FmCopy(fm1, fromf ? f0x1 : fax1);
  if (trans_on) {
    FmCopy(fm2, fromf ? f0x2 : fax2);
    if (trans_x3) {FmCopy(fm3, fromf ? f0x3 : fax3);}
  }
  if (t2km1.extent_int(0) > 0 && t2k1.extent_int(1) >= M1_T2_F1 + 3) {
    auto k1_ = t2k1;
    auto km_ = t2km1;
    par_for("m1_fm_k1init", DevExeSpace(), 0, t2km1.extent_int(0) - 1, 0, 2,
            0, t2km1.extent_int(2) - 1, 0, t2km1.extent_int(3) - 1,
            0, t2km1.extent_int(4) - 1,
    KOKKOS_LAMBDA(const int m, const int d, const int k, const int j, const int i) {
      km_(m,d,k,j,i) = k1_(m,M1_T2_F1+d,k,j,i);
    });
  }
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1> implicit_hr_cmem_fm: F_m set from the "
              << (fromf ? "face flux f0x*" : "applied face fluxes fax*") << " at cycle "
              << pmy_pack->pmesh->ncycle << std::endl;
  }
}

//! the hesdirk2 stage copies (as t2f* := f0x* and back)
void RadiationM1::FmSave() {
  if (!impl_cm_fm) return;
  FmEnsure();
  FmCopy(t2fm1, fm1);
  if (trans_on) {
    FmCopy(t2fm2, fm2);
    if (trans_x3) {FmCopy(t2fm3, fm3);}
  }
}

void RadiationM1::FmLoad() {
  if (!impl_cm_fm) return;
  FmEnsure();
  FmCopy(fm1, t2fm1);
  if (trans_on) {
    FmCopy(fm2, t2fm2);
    if (trans_x3) {FmCopy(fm3, t2fm3);}
  }
}

//! the stage-2 start: the Heun (or w0/w1) average of U^n and Y2, as avgf of f0x*
void RadiationM1::FmAvg(bool heun, Real w0, Real w1) {
  if (!impl_cm_fm) return;
  FmEnsure();
  auto avg = [&](DvceArray4D<Real> &a, DvceArray4D<Real> &b) {
    if (a.extent_int(0) == 0 || b.extent(0) != a.extent(0)) return;
    auto a_ = a;
    auto b_ = b;
    par_for("m1_fm_avg", DevExeSpace(), 0, a.extent_int(0) - 1, 0, a.extent_int(1) - 1,
            0, a.extent_int(2) - 1, 0, a.extent_int(3) - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      a_(m,k,j,i) = heun ? 0.5*(b_(m,k,j,i) + a_(m,k,j,i))
                         : (w1*b_(m,k,j,i) + w0*a_(m,k,j,i));
    });
  };
  avg(fm1, t2fm1);
  if (trans_on) {
    avg(fm2, t2fm2);
    if (trans_x3) {avg(fm3, t2fm3);}
  }
}

//! the stage increment, as t2inc's F channels (m1_t2_inc + the posfix take-back pcomp;
//! time2_face_fct does not scale it: F_m is no E transport, and the cell-F increment
//! csinc it must stay consistent with is not scaled either)
void RadiationM1::FmFormInc(bool two, Real c1, Real c2, bool pcomp, Real ccomp) {
  if (!impl_cm_fm) return;
  FmEnsure();
  if (t2incm.extent_int(0) == 0) return;
  auto k1_ = t2km1;
  auto k2_ = t2km2;
  auto in_ = t2incm;
  const int nmb1 = t2incm.extent_int(0) - 1;
  const int n3 = t2incm.extent_int(2), n2 = t2incm.extent_int(3);
  const int n1 = t2incm.extent_int(4);
  par_for("m1_fm_inc", DevExeSpace(), 0, nmb1, 0, 2, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int d, const int k, const int j, const int i) {
    Real v = c1*k1_(m,d,k,j,i);
    if (two) {v += c2*k2_(m,d,k,j,i);}
    in_(m,d,k,j,i) = v;
  });
  if (pcomp && t2pm.extent_int(0) > nmb1) {
    auto pm_ = t2pm;
    const bool d2 = (n2 > 1), d3 = (n3 > 1);
    par_for("m1_fm_pcomp", DevExeSpace(), 0, nmb1, 0, 2, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int d, const int k, const int j, const int i) {
      Real w = pm_(m,k,j,i);
      if (d == 0 && i > 0) {w = fmax(w, pm_(m,k,j,i-1));}
      if (d == 1) {w = d2 ? ((j > 0) ? fmax(w, pm_(m,k,j-1,i)) : w) : 0.0;}
      if (d == 2) {w = d3 ? ((k > 0) ? fmax(w, pm_(m,k-1,j,i)) : w) : 0.0;}
      if (w > 0.0) {in_(m,d,k,j,i) += ccomp*k1_(m,d,k,j,i);}
    });
  }
}

//! Time2RotCorrect: increment = c K1
void RadiationM1::FmIncFromK1(Real c) {
  if (!impl_cm_fm) return;
  FmEnsure();
  if (t2incm.extent_int(0) == 0) return;
  auto k1_ = t2km1;
  auto in_ = t2incm;
  par_for("m1_fm_rcinc", DevExeSpace(), 0, t2incm.extent_int(0) - 1, 0, 2,
          0, t2incm.extent_int(2) - 1, 0, t2incm.extent_int(3) - 1,
          0, t2incm.extent_int(4) - 1,
  KOKKOS_LAMBDA(const int m, const int d, const int k, const int j, const int i) {
    in_(m,d,k,j,i) = c*k1_(m,d,k,j,i);
  });
}

//! Time2PosFix: drop the stage-1 increment of the faces with a marked side
void RadiationM1::FmPfixDrop() {
  if (!impl_cm_fm || t2incm.extent_int(0) == 0 || t2pm.extent_int(0) == 0) return;
  auto in_ = t2incm;
  auto pm_ = t2pm;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const bool d2 = (indcs.nx2 > 1), d3 = (indcs.nx3 > 1);
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n3 = t2incm.extent_int(2), n2 = t2incm.extent_int(3);
  const int n1 = t2incm.extent_int(4);
  par_for("m1_fm_pdrop", DevExeSpace(), 0, nmb1, 0, 2, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int d, const int k, const int j, const int i) {
    Real w = pm_(m,k,j,i);
    if (d == 0 && i > 0) {w = fmax(w, pm_(m,k,j,i-1));}
    if (d == 1) {w = d2 ? ((j > 0) ? fmax(w, pm_(m,k,j-1,i)) : w) : 0.0;}
    if (d == 2) {w = d3 ? ((k > 0) ? fmax(w, pm_(m,k-1,j,i)) : w) : 0.0;}
    if (w > 0.0) {in_(m,d,k,j,i) = 0.0;}
  });
}

void RadiationM1::FmK1FromK2() {
  if (!impl_cm_fm || t2km1.extent_int(0) == 0) return;
  Kokkos::deep_copy(DevExeSpace(), t2km1, t2km2);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::FmStepStart
//! \brief the solve start, as f0x*: the cubed-sphere seam average of the transverse
//! faces (CubedSeamFaceAverage), the old copy fm*n, and under a hesdirk2 stage solve the
//! stage increment (as m1_t2_f1n / f2n / f3n)

void RadiationM1::FmStepStart(bool t2st) {
  if (!impl_cm_fm) return;
  FmEnsure();
  if (cs_geom && trans_on && trans_x3) {CubedSeamFaceAverageOf(fm2, fm3, false);}
  Kokkos::deep_copy(DevExeSpace(), fm1n, fm1);
  if (trans_on) {
    Kokkos::deep_copy(DevExeSpace(), fm2n, fm2);
    if (trans_x3) {Kokkos::deep_copy(DevExeSpace(), fm3n, fm3);}
  }
  if (!t2st || t2incm.extent_int(0) == 0) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto in_ = t2incm;
  auto f1n_ = fm1n;
  par_for("m1_fm_f1n", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    f1n_(m,k,j,i) += in_(m,0,k,j,i);
  });
  if (trans_on) {
    auto f2n_ = fm2n;
    par_for("m1_fm_f2n", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      f2n_(m,k,j,i) += in_(m,1,k,j,i);
    });
    if (trans_x3) {
      auto f3n_ = fm3n;
      par_for("m1_fm_f3n", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        f3n_(m,k,j,i) += in_(m,2,k,j,i);
      });
    }
  }
}

//! the face slope of the solve (stage 1 -> K2, otherwise K1), as m1_t2_kf1 / 2 / 3
void RadiationM1::FmSlope(bool stage1, Real dt) {
  if (!impl_cm_fm || t2km1.extent_int(0) == 0) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto kk_ = stage1 ? t2km2 : t2km1;
  const Real fk = 1.0/dt;
  auto f1_ = fm1, f1n_ = fm1n;
  par_for("m1_fm_kf1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    kk_(m,0,k,j,i) = (f1_(m,k,j,i) - f1n_(m,k,j,i))*fk;
  });
  if (trans_on) {
    auto f2_ = fm2, f2n_ = fm2n;
    par_for("m1_fm_kf2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      kk_(m,1,k,j,i) = (f2_(m,k,j,i) - f2n_(m,k,j,i))*fk;
    });
    if (trans_x3) {
      auto f3_ = fm3, f3n_ = fm3n;
      par_for("m1_fm_kf3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        kk_(m,2,k,j,i) = (f3_(m,k,j,i) - f3n_(m,k,j,i))*fk;
      });
    }
  }
}

}  // namespace radm1
