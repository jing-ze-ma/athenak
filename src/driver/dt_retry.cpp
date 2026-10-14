//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file dt_retry.cpp
//! \brief fin4-dtretry: the cycle checkpoint and the dt cut-and-retry (dt_retry.hpp)

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "hydro/hydro.hpp"
#include "rad_m1/rad_m1.hpp"
#include "hydro/hydro_dtretry_visit.hpp"
#include "rad_m1/rad_m1_dtretry_visit.hpp"
#include "driver/dt_retry.hpp"

namespace dtretry {

namespace {

// the evolving arrays of dt_retry_views = minimal: what a cycle reads from the cycle
// before (the restart file's set, plus the lagged closure/opacity inputs); every other
// View is scratch that a cycle rebuilds before reading it (the audit, dt_retry_audit,
// lists the Views a failed attempt changed; dt_retry_poison checks the rest are scratch)
const char *kDtrMinimal[] = {
  // Hydro
  "hyd.u0", "hyd.w0", "hyd.wtemp", "hyd.wder", "hyd.eint_rst", "hyd.c2p_acc",
  "hyd.c2p_ref_e", "hyd.c2p_ref_d",
  // RadiationM1
  "m1.u0", "m1.f0x1", "m1.f0x2", "m1.f0x3", "m1.ipred", "m1.ipred2", "m1.t2k1",
  "m1.vet_prev", "m1.vsc_d0", "m1.vsc_d1", "m1.ctr_mem", "m1.csk1", "m1.csk2",
  "m1.csinc", "m1.vra_d", "m1.opac", "m1.tau_ten", "m1.vet_cell", "m1.vgd_rst",
  "m1.vlat_cs",
};

// copy-assign anything, arrays element by element
template <class T> void Assign(T &a, const T &b) {
  if constexpr (std::is_array_v<T>) {
    for (std::size_t q = 0; q < std::extent_v<T>; ++q) {Assign(a[q], b[q]);}
  } else {
    a = b;
  }
}
template <class T> struct Holder {T x;};

template <class T>
struct ValEntry : public Entry {
  T *p;
  Holder<T> c;
  explicit ValEntry(T *pp) : p(pp) {}
  void Save() override {Assign(c.x, *p);}
  void Restore() override {Assign(*p, c.x);}
};

template <class T>
struct ViewEntry : public Entry {
  using Clone = Kokkos::View<typename T::non_const_data_type, typename T::array_layout,
                             typename T::memory_space>;
  T *p;
  T orig;
  Clone clone;
  bool deep, have = false;
  ViewEntry(T *pp, bool d) : p(pp), deep(d) {is_view = true;}
  bool Same(const T &a, const Clone &b) const {
    if (a.data() == nullptr || b.data() == nullptr) return false;
    for (int r = 0; r < static_cast<int>(T::rank); ++r) {
      if (a.extent(r) != b.extent(r)) return false;
    }
    return true;
  }
  void Save() override {
    orig = *p;
    have = false;
    if (!deep || p->data() == nullptr || p->span() == 0) return;
    if (!Same(*p, clone)) {
      clone = Clone(Kokkos::view_alloc(std::string("dtr:") + name,
                                       Kokkos::WithoutInitializing), p->layout());
    }
    Kokkos::deep_copy(clone, *p);
    have = true;
  }
  void Restore() override {
    bool moved = (p->data() != orig.data());
    for (int r = 0; r < static_cast<int>(T::rank) && !moved; ++r) {
      moved = (p->extent(r) != orig.extent(r));
    }
    if (moved) {*p = orig;}
    if (have) {Kokkos::deep_copy(*p, clone);}
  }
  std::size_t Bytes() const override {
    return have ? clone.span()*sizeof(typename T::non_const_value_type) : 0;
  }
  bool Differs() override {
    if (!have) return false;
    if (p->data() != orig.data()) return true;
    auto a = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), *p);
    auto b = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), clone);
    return std::memcmp(a.data(), b.data(),
                       a.span()*sizeof(typename T::non_const_value_type)) != 0;
  }
  void Poison() override {
    if constexpr (std::is_floating_point_v<typename T::non_const_value_type>) {
      if (p->data() != nullptr && p->span() > 0) {
        Kokkos::deep_copy(*p, std::numeric_limits<
                              typename T::non_const_value_type>::quiet_NaN());
      }
    }
  }
};

template <class T>
struct DualEntry : public Entry {
  T *p;
  T orig;
  ViewEntry<typename T::t_dev> d;
  ViewEntry<typename T::t_host> h;
  bool nsh = false, nsd = false;
  DualEntry(T *pp, bool deep) : p(pp), d(&pp->d_view, deep), h(&pp->h_view, deep) {
    is_view = true;
  }
  void Save() override {
    orig = *p;
    d.p = &p->d_view; h.p = &p->h_view;
    d.name = name + ".d"; h.name = name + ".h";
    d.Save(); h.Save();
    nsh = p->need_sync_host(); nsd = p->need_sync_device();
  }
  void Restore() override {
    if (p->d_view.data() != orig.d_view.data() || p->h_view.data() != orig.h_view.data()) {
      *p = orig;
    }
    d.p = &p->d_view; h.p = &p->h_view;
    d.Restore(); h.Restore();
    p->clear_sync_state();
    if (nsd) {p->modify_host();}
    if (nsh) {p->modify_device();}
  }
  std::size_t Bytes() const override {return d.Bytes() + h.Bytes();}
};

// the visitor that registers every member of one module (prefix "hyd." / "m1.")
struct Binder {
  std::vector<std::unique_ptr<Entry>> *ent;
  std::string pre;
  const std::set<std::string> *sel;
  bool all;
  bool Deep(const std::string &n) const {
    return all || sel->count(pre + n) > 0;
  }
  template <class T> void value(const char *n, T &x) {
    auto e = std::make_unique<ValEntry<T>>(&x);
    e->name = pre + n;
    ent->push_back(std::move(e));
  }
  template <class T> void view(const char *n, T &x) {
    auto e = std::make_unique<ViewEntry<T>>(&x, Deep(n));
    e->name = pre + n;
    ent->push_back(std::move(e));
  }
  template <class T> void dual(const char *n, T &x) {
    auto e = std::make_unique<DualEntry<T>>(&x, Deep(n));
    e->name = pre + n;
    ent->push_back(std::move(e));
  }
};

}  // namespace

//----------------------------------------------------------------------------------------
DtRetry::DtRetry(ParameterInput *pin, Mesh *pm) {
  const bool fresh = !global_variable::restart_run;
  if (pin->DoesParameterExist("time", "dt_retry")) {
    on = pin->GetBoolean("time", "dt_retry");
  } else {
    on = pin->GetOrAddBoolean("time", "dt_retry", fresh);
  }
  if (!on) return;
  // the wiring covers Hydro + RadiationM1 (the implicit scheme); anything else: off
  auto *pk = pm->pmb_pack;
  if (pk->pmhd != nullptr || pk->prad != nullptr || pk->pz4c != nullptr ||
      pk->padm != nullptr || pk->ppart != nullptr || pk->pturb != nullptr ||
      pk->pradm1 == nullptr || pm->adaptive || pm->multilevel || pk->pradm1->vgd_async) {
    if (global_variable::my_rank == 0) {
      std::cout << "<time> dt_retry: not wired for this module set (wired: <hydro> + "
                << "<rad_m1>, no mesh refinement, no vet_gd_async) -> off" << std::endl;
    }
    on = false;
    pin->SetBoolean("time", "dt_retry", false);
    return;
  }
  fac_ = pin->GetOrAddReal("time", "dt_retry_fac", 0.5);
  grow_ = pin->GetOrAddReal("time", "dt_retry_grow", 1.25);
  maxcut_ = pin->GetOrAddInteger("time", "dt_retry_max", 8);
  ngood_req_ = pin->GetOrAddInteger("time", "dt_retry_ngood", 10);
  std::string vm = pin->GetOrAddString("time", "dt_retry_views", "minimal");
  all_views_ = (vm == "all");
  if (!all_views_ && vm != "minimal") {
    std::cout << "### FATAL ERROR: <time>/dt_retry_views must be minimal or all" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (pin->DoesParameterExist("time", "dt_retry_views_add")) {
    std::stringstream ss(pin->GetString("time", "dt_retry_views_add"));
    std::string w;
    while (std::getline(ss, w, ',')) {if (!w.empty()) extra_views_.push_back(w);}
  }
  // debug / audit keys, read only when named (the parameter dump is unchanged otherwise)
  if (pin->DoesParameterExist("time", "dt_retry_force_cycle")) {
    force_cycle_ = pin->GetInteger("time", "dt_retry_force_cycle");
    force_fac_ = pin->GetOrAddReal("time", "dt_retry_force_fac", 0.5);
  }
  if (pin->DoesParameterExist("time", "dt_retry_audit")) {
    audit_ = pin->GetBoolean("time", "dt_retry_audit");
  }
  if (pin->DoesParameterExist("time", "dt_retry_poison")) {
    poison_ = pin->GetInteger("time", "dt_retry_poison");
    if (poison_ > 0 && !all_views_) {
      std::cout << "### FATAL ERROR: <time>/dt_retry_poison needs dt_retry_views = all"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // the cap and the good-cycle count of the straight run (a restart file's dump)
  if (pin->DoesParameterExist("time", "dt_retry_cap")) {
    pm->dt_cap = pin->GetReal("time", "dt_retry_cap");
  }
  if (pin->DoesParameterExist("time", "dt_retry_ngood_now")) {
    ngood_ = pin->GetInteger("time", "dt_retry_ngood_now");
  }
  pk->pradm1->dtr_on = true;
  if (global_variable::my_rank == 0) {
    std::cout << "<time> dt_retry on: fac " << fac_ << " max " << maxcut_ << " ngood "
              << ngood_req_ << " grow " << grow_ << " views "
              << (all_views_ ? "all" : "minimal") << " cap "
              << pm->dt_cap << std::endl;
  }
}

//----------------------------------------------------------------------------------------
void DtRetry::Build(Mesh *pm) {
  std::set<std::string> &sel = sel_;
  for (const char *s : kDtrMinimal) {sel.insert(s);}
  for (auto &s : extra_views_) {sel.insert(s);}
  auto *pk = pm->pmb_pack;
  if (pk->phydro != nullptr) {
    Binder b{&ent_, "hyd.", &sel, all_views_};
    pk->phydro->DtRetryVisit(b);
  }
  {
    Binder b{&ent_, "m1.", &sel, all_views_};
    pk->pradm1->DtRetryVisit(b);
  }
  {
    auto e = std::make_unique<ValEntry<EventCounters>>(&(pm->ecounter));
    e->name = "mesh.ecounter";
    ent_.push_back(std::move(e));
  }
  built_ = true;
}

//----------------------------------------------------------------------------------------
void DtRetry::Save(Mesh *pm) {
  if (!on) return;
  if (!built_) {Build(pm);}
  Kokkos::Timer t;
  pm->pmb_pack->pradm1->DtRetryQuiesce(false);
  for (auto &e : ent_) {e->Save();}
  Kokkos::fence();
  tsave_ += t.seconds();
  if (nsave_ == 0) {
    std::size_t b = 0;
    int nv = 0, nd = 0;
    for (auto &e : ent_) {
      b += e->Bytes();
      if (e->is_view) {nv += 1; if (e->Bytes() > 0) nd += 1;}
    }
    bytes_dev_ = b;
    if (global_variable::my_rank == 0) {
      std::cout << "<time> dt_retry checkpoint: " << ent_.size() << " members, " << nv
                << " Views (" << nd << " copied, " << static_cast<double>(b)/1.0e6
                << " MB on rank 0):";
      for (auto &e : ent_) {
        if (e->Bytes() > 0) {
          std::cout << " " << e->name << " " << static_cast<double>(e->Bytes())/1.0e6;
        }
      }
      std::cout << std::endl;
    }
  }
  nsave_ += 1;
}

//----------------------------------------------------------------------------------------
bool DtRetry::Failed(Mesh *pm) {
  if (!on) return false;
  auto *m1 = pm->pmb_pack->pradm1;
  if (force_cycle_ >= 0 && !forced_done_ && pm->ncycle == force_cycle_) {
    forced_done_ = true;
    m1->dtr_fail = true;
    m1->dtr_why = 9;
  }
  return m1->dtr_fail;
}

//----------------------------------------------------------------------------------------
bool DtRetry::Retry(Mesh *pm, Driver *pd) {
  auto *m1 = pm->pmb_pack->pradm1;
  const int why = m1->dtr_why;
  const Real res = m1->dtr_resid;
  const Real dto = pm->dt;
  const Real f = (why == 9) ? force_fac_ : fac_;
  ncut_ += 1;
  nretry_ += 1;
  if (ncut_ > maxcut_) {
    if (global_variable::my_rank == 0) {
      std::cout << "### FATAL ERROR: <time> dt_retry: cycle " << pm->ncycle << " failed "
                << ncut_ << " times (dt_retry_max " << maxcut_ << "), dt " << dto
                << std::endl;
    }
    return false;
  }
  // audit: which Views did the failed attempt change (before they are put back)
  std::vector<Entry*> changed;
  if (audit_ || poison_ > 0) {
    for (auto &e : ent_) {
      if (e->is_view && e->Differs()) {changed.push_back(e.get());}
    }
  }
  m1->DtRetryQuiesce(true);
  for (auto &e : ent_) {e->Restore();}
  Kokkos::fence();
  m1->DtRetryAfterRestore();
  m1->dtr_fail = false;
  m1->dtr_why = 0;
  m1->dtr_resid = 0.0;
  if (poison_ > 0 && why == 9) {
    // dt_retry_views = all: every View the failed attempt changed that the minimal list
    // does NOT checkpoint goes to NaN after the restore; the retried cycle is bitwise the
    // straight one only if each of them is scratch (written before it is read)
    int np = 0;
    for (auto *e : changed) {
      if (sel_.count(e->name) == 0) {e->Poison(); np += 1;}
    }
    Kokkos::fence();
    if (global_variable::my_rank == 0) {
      std::cout << "<time> dt_retry poison: " << np << " changed Views outside the minimal "
                << "list set to NaN" << std::endl;
    }
  }
  if (audit_ && global_variable::my_rank == 0) {
    std::cout << "<time> dt_retry audit cycle " << pm->ncycle << ": " << changed.size()
              << " checkpointed Views changed by the failed attempt:";
    for (auto *e : changed) {std::cout << " " << e->name;}
    std::cout << std::endl;
  }
  Real dtn = f*dto;
  if (why != 9 || f != 1.0) {
    pm->dt_cap = dtn;
    ngood_ = 0;
  }
  if (dtn != dto && pm->pmb_pack->pradm1 != nullptr) {
    dtn = pm->pmb_pack->pradm1->VetGdAlignDt(pm->time, dtn);
  }
  pm->dt = dtn;
  static const char *whys[] = {"?", "NON-CONVERGED", "DIVERGED", "DIVERGED(masked)",
                               "eint<=0", "?", "?", "?", "?", "forced"};
  if (global_variable::my_rank == 0) {
    std::cout << "DT_RETRY cycle=" << pm->ncycle << " t=" << pm->time << " dt " << dto
              << " -> " << dtn << " reason=" << whys[std::min(std::max(why, 0), 9)]
              << " resid=" << res << " cut " << ncut_ << std::endl;
  }
  (void) pd;
  return true;
}

//----------------------------------------------------------------------------------------
void DtRetry::Accept(Mesh *pm) {
  if (!on) return;
  ncut_ = 0;
  if (pm->dt_cap > 0.0) {
    ngood_ += 1;
    if (ngood_ > ngood_req_) {pm->dt_cap *= grow_;}
  }
}

//----------------------------------------------------------------------------------------
void DtRetry::AfterNewDt(Mesh *pm, ParameterInput *pin) {
  if (!on) return;
  if (pm->dt_cap > 0.0 && pm->dt_cap >= pm->dt_uncapped) {
    if (global_variable::my_rank == 0) {
      std::cout << "DT_RETRY recovered cycle=" << pm->ncycle << " t=" << pm->time
                << " cap " << pm->dt_cap << " >= CFL dt " << pm->dt_uncapped
                << ": cap off" << std::endl;
    }
    pm->dt_cap = -1.0;
    ngood_ = 0;
  }
  // recorded for a restart file (written after this) only while a cap is or was set
  if (pm->dt_cap > 0.0 || pin->DoesParameterExist("time", "dt_retry_cap")) {
    pin->SetReal("time", "dt_retry_cap", pm->dt_cap);
    pin->SetInteger("time", "dt_retry_ngood_now", ngood_);
  }
}

//----------------------------------------------------------------------------------------
void DtRetry::Report(Mesh *pm) {
  if (!on || global_variable::my_rank != 0) return;
  std::cout << "<time> dt_retry: " << nsave_ << " checkpoints (" << tsave_ << " s, "
            << static_cast<double>(bytes_dev_)/1.0e6 << " MB on rank 0), " << nretry_
            << " retries, cap " << pm->dt_cap << std::endl;
}

}  // namespace dtretry
