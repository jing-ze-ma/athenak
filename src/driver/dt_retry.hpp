#ifndef DRIVER_DT_RETRY_HPP_
#define DRIVER_DT_RETRY_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file dt_retry.hpp
//! \brief fin4-dtretry: reject a failed cycle and take it again with a smaller dt.
//!
//! <time>/dt_retry (DEFAULT true on a fresh run; a restart whose input lacks the key keeps
//! false, as time2_enth_be): at the top of every cycle the Driver checkpoints the evolving
//! state of the Hydro and the RadiationM1 modules (DtRetry::Save).  A cycle in which a
//! module raises a failure (RadiationM1::dtr_fail: a Picard solve that ends DIVERGED, or
//! NON-CONVERGED outside a hesdirk2 stage, whose own redo is the backward-Euler step) is
//! put back exactly (DtRetry::Restore) and taken again with dt x dt_retry_fac; the cut dt
//! becomes a cap that, after dt_retry_ngood accepted cycles in a row, grows by
//! dt_retry_grow per cycle and is dropped once the CFL dt is below it.  The cap and the
//! good-cycle count are kept in <time>/dt_retry_cap, dt_retry_ngood_now (the parameter
//! dump of a restart file), so a restart takes the step the straight run took.
//! A run that never triggers it is bitwise the run without it: the checkpoint only copies.
//!
//! The checkpoint is the member visitor generated from clang's AST (tools/dtretry_gen.py:
//! rad_m1/rad_m1_dtretry_visit.hpp, hydro/hydro_dtretry_visit.hpp): every non-View member
//! is copied by value, and every View member in the selected set is deep-copied into a
//! clone of the same layout and memory space (its handle is kept too, so a View swapped
//! or reallocated during the failed attempt is put back).  dt_retry_views = minimal
//! (default: the evolving arrays, kDtrMinimal in dt_retry.cpp) | all (every View: memory
//! x ~2, a debug and audit mode).

#include <memory>
#include <set>
#include <string>
#include <vector>

#include "athena.hpp"

class Mesh;
class Driver;
class ParameterInput;

namespace dtretry {

// one checkpointed member
struct Entry {
  std::string name;
  bool is_view = false;
  virtual ~Entry() = default;
  virtual void Save() = 0;
  virtual void Restore() = 0;
  virtual std::size_t Bytes() const {return 0;}
  virtual bool Differs() {return false;}   // audit: changed since Save (Views only)
  virtual void Poison() {}                 // audit: fill a View with NaN
};

// fin51-r-fix: file-scope state of a problem generator that a cycle reads from the cycle
// before (or that a once-per-cycle update guards with a cycle stamp) joins the checkpoint:
// the pgen registers it once at setup (RegisterValue / RegisterView, below) and
// DtRetry::Build adds it as "pgen.<name>" entries (Views always deep-copied).  Without it
// a retried cycle keeps what the failed attempt left there (he_star_m1: the frozen-MLT
// ramp weight hs_mw_ was already at w(t) on the retry, so esrc kept w(t) where the straight
// run's start-of-cycle solves still saw w(t - dt)).
std::vector<std::unique_ptr<Entry>> &PgenEntries();

template <class T>
struct PgenValEntry : public Entry {
  T *p;
  T c{};
  explicit PgenValEntry(T *pp) : p(pp) {}
  void Save() override {c = *p;}
  void Restore() override {*p = c;}
};

template <class V>
struct PgenViewEntry : public Entry {
  using Clone = Kokkos::View<typename V::non_const_data_type, typename V::array_layout,
                             typename V::memory_space>;
  V *p;
  V orig;
  Clone clone;
  bool have = false;
  explicit PgenViewEntry(V *pp) : p(pp) {is_view = true;}
  void Save() override {
    orig = *p;
    have = false;
    if (p->data() == nullptr || p->span() == 0) return;
    bool same = (clone.data() != nullptr);
    for (int r = 0; r < static_cast<int>(V::rank) && same; ++r) {
      same = (clone.extent(r) == p->extent(r));
    }
    if (!same) {
      clone = Clone(Kokkos::view_alloc(std::string("dtr:") + name,
                                       Kokkos::WithoutInitializing), p->layout());
    }
    Kokkos::deep_copy(clone, *p);
    have = true;
  }
  void Restore() override {
    if (p->data() != orig.data()) {*p = orig;}
    if (have) {Kokkos::deep_copy(*p, clone);}
  }
  std::size_t Bytes() const override {
    return have ? clone.span()*sizeof(typename V::non_const_value_type) : 0;
  }
};

template <class T>
void RegisterValue(const char *name, T *p) {
  auto e = std::make_unique<PgenValEntry<T>>(p);
  e->name = std::string("pgen.") + name;
  PgenEntries().push_back(std::move(e));
}
template <class V>
void RegisterView(const char *name, V *p) {
  auto e = std::make_unique<PgenViewEntry<V>>(p);
  e->name = std::string("pgen.") + name;
  PgenEntries().push_back(std::move(e));
}

class DtRetry {
 public:
  DtRetry(ParameterInput *pin, Mesh *pm);
  bool on = false;
  // checkpoint at the top of a cycle (no-op when off)
  void Save(Mesh *pm);
  // a module raised a failure in this cycle (MPI-consistent: the flags are set after
  // reductions on every rank)
  bool Failed(Mesh *pm);
  // put the checkpoint back and cut dt; false = no further cut allowed (caller stops)
  bool Retry(Mesh *pm, Driver *pd);
  // an accepted cycle, before Mesh::NewTimeStep: the cap grows after ngood cycles
  void Accept(Mesh *pm);
  // after Mesh::NewTimeStep: the cap is dropped once the CFL dt is below it; recorded
  void AfterNewDt(Mesh *pm, ParameterInput *pin);
  // end of run: checkpoints, retries, checkpoint time
  void Report(Mesh *pm);

 private:
  void Build(Mesh *pm);
  bool built_ = false;
  bool all_views_ = false;
  bool audit_ = false;
  int poison_ = 0;
  Real fac_ = 0.5, grow_ = 1.25;
  int maxcut_ = 8, ngood_req_ = 10;
  int ngood_ = 0;        // accepted cycles since the last cut
  int ncut_ = 0;         // cuts of the current cycle
  int force_cycle_ = -1; // debug: raise a failure at this cycle (once)
  Real force_fac_ = 0.5;
  bool forced_done_ = false;
  std::vector<std::string> extra_views_;
  std::set<std::string> sel_;
  std::vector<std::string> poison_match_;
  std::vector<std::unique_ptr<Entry>> ent_;
  std::size_t bytes_dev_ = 0;
  double tsave_ = 0.0;
  long nsave_ = 0, nretry_ = 0;
};

}  // namespace dtretry
#endif  // DRIVER_DT_RETRY_HPP_
