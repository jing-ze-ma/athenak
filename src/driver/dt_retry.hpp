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

#include <cstdint>
#include <memory>
#include <set>
#include <string>
#include <utility>
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
  virtual bool OnHost() const {return false;}   // gpumem-1022: copy in pinned host memory
  virtual bool Differs() {return false;}   // audit: changed since Save (Views only)
  virtual void Poison() {}                 // audit: fill a View with NaN
};

// fin51-r-fix: file-scope state of a problem generator that a cycle reads from the cycle
// before (or that a once-per-cycle update guards with a cycle stamp) joins the
// checkpoint: the pgen registers it once at setup (RegisterValue / RegisterView, below)
// and DtRetry::Build adds it as "pgen.<name>" entries (Views always deep-copied).
// Without it a retried cycle keeps what the failed attempt left there (he_star_m1: the
// frozen-MLT ramp weight hs_mw_ was already at w(t) on the retry, so esrc kept w(t)
// where the straight run's start-of-cycle solves still saw w(t - dt)).
std::vector<std::unique_ptr<Entry>> &PgenEntries();
// pgrst-1021: a second, independent copy of every registered pgen member (the dt_retry
// checkpoint above takes the PgenEntries list over at its Build).  The driver saves it
// at the start of the hesdirk2 stage loop and puts it back with the backward-Euler redo
// (Time2Restore): the redo reruns the hydro stages from U^n, so the pgen registers that
// advance in every stage (ry_per_accretor acc0_/acc1_/stage_ctr_, tracer + L_z
// registers) must restart from their stage-loop values too, else the failed attempt is
// booked twice and the stage counter goes out of phase.  Nothing a kernel reads in a
// straight (no-redo) cycle is touched: the default path is unchanged.
std::vector<std::unique_ptr<Entry>> &PgenStageEntries();
void PgenStageSave();
void PgenStageRestore();
void PgenStageClear();   // the Driver drops the copies before Kokkos::finalize
// pgrst-posfix-1021: a third copy, taken before the stage-1 radiation solve when
// <rad_m1>/time2_posfix is on and put back before each local re-solve (the M1 post hook
// of the discarded attempt booked into the registers)
std::vector<std::unique_ptr<Entry>> &PgenSolveEntries();
void PgenSolveSave();
void PgenSolveRestore();

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
  auto f = std::make_unique<PgenValEntry<T>>(p);   // pgrst-1021: stage-loop copy
  f->name = std::string("pgen.") + name;
  PgenStageEntries().push_back(std::move(f));
  auto g = std::make_unique<PgenValEntry<T>>(p);   // pgrst-posfix-1021: solve copy
  g->name = std::string("pgen.") + name;
  PgenSolveEntries().push_back(std::move(g));
}
template <class V>
void RegisterView(const char *name, V *p) {
  auto e = std::make_unique<PgenViewEntry<V>>(p);
  e->name = std::string("pgen.") + name;
  PgenEntries().push_back(std::move(e));
  auto f = std::make_unique<PgenViewEntry<V>>(p);  // pgrst-1021: stage-loop copy
  f->name = std::string("pgen.") + name;
  PgenStageEntries().push_back(std::move(f));
  auto g = std::make_unique<PgenViewEntry<V>>(p);  // pgrst-posfix-1021: solve copy
  g->name = std::string("pgen.") + name;
  PgenSolveEntries().push_back(std::move(g));
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
  bool host_ = false;      // gpumem-1022 <time>/dt_retry_host: clones in pinned host memory
  int poison_ = 0;
  Real fac_ = 0.5, grow_ = 1.25;
  int maxcut_ = 8, ngood_req_ = 10;
  int ngood_ = 0;        // accepted cycles since the last cut
  int ncut_ = 0;         // cuts of the current cycle
  int force_cycle_ = -1; // debug: raise a failure at this cycle (once)
  Real force_fac_ = 0.5;
  int t2_be_ = 2;          // t2befb-1020: <time>/dt_retry_t2_be (time2_reject cuts per
                           // cycle before the backward-Euler redo; < 0 = never)
  int nt2cut_ = 0;         // t2befb-1020: time2_reject cuts of the current cycle
  long nt2be_ = 0;         // t2befb-1020: cycles that took the backward-Euler redo
  bool t2streak_ = false;  // t2befb-1020: the last cycle took it (the next one is armed)
  // prod-retry-1022: this cycle already went to the backward-Euler redo (t2fb): the
  // restore of a later failed attempt must not disarm it (it resets dtr_t2be)
  bool t2be_cyc_ = false;
  Real dt0_ = 0.0, cap0_ = -1.0;  // t2befb-1020: the cycle's dt and cap at Save
  int ngood0_ = 0;
  Real t2_fac_ = -1.0;     // t2retry-1020: <time>/dt_retry_fac_t2 (reason time2_reject; < 0 = fac_)
  bool forced_done_ = false;
  std::vector<std::string> extra_views_;
  std::set<std::string> sel_;
  std::vector<std::string> poison_match_;
  std::vector<std::unique_ptr<Entry>> ent_;
  std::size_t bytes_dev_ = 0;
  std::size_t bytes_host_ = 0;   // gpumem-1022: of bytes_dev_, the part on the host
  double tsave_ = 0.0;
  long nsave_ = 0, nretry_ = 0;
};

// m1ctr-rollback-1017: the EVENT COUNTERS of an attempt the code throws away are thrown
// away with it, so the event log (and the m1 positivity totals) report what the state
// kept.  Discarded attempts: a dt_retry cycle (DtRetry::Retry puts the state back, but
// not Mesh::ecounter nor the <rad_m1> device counters pos_cnt_d / acl_cnt_d /
// opn_nskip_d), a hesdirk2 step redone with backward Euler (Time2Restore: U^n back, the
// hydro stages run again), and the hesdirk2 helper solves whose result is put back
// (Time2StartProbe, Time2RotCorrect).  Counters only: nothing a kernel reads is touched,
// so the evolution is bitwise the same with it on or off.  <time>/counters_rollback
// (default true; read only when named, the parameter dump is unchanged).  Slot 0 = the
// cycle checkpoint (dt_retry), slot 1 = the stage start (backward-Euler redo).
class CounterKeep {
 public:
  CounterKeep(ParameterInput *pin, Mesh *pm);
  bool on = true;
  void Save(Mesh *pm, int slot);
  void Restore(Mesh *pm, int slot, int why);   // why: 0 dt_retry, 1 BE redo, 2 helper
  void Report(Mesh *pm);

 private:
  struct Snap {
    bool have = false;
    Real ec[3] = {0, 0, 0};               // Mesh::ecounter efloor_de, vceil_de, dfloor_ke
    std::int64_t ni[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};   // its integer counters
    DvceArray1D<Real> pos, acl, opn;
    bool hpos = false, hacl = false, hopn = false;
  };
  Snap s_[3];   // slot 2: t2pos-1020 the stage-1 start of a time2_posfix re-solve
  long nrb_[3] = {0, 0, 0};
  double dun_[3] = {0.0, 0.0, 0.0};    // E-floor created energy discarded (m1 pos UN)
  double dvc_[3] = {0.0, 0.0, 0.0};    // vceil_de discarded
};

}  // namespace dtretry
#endif  // DRIVER_DT_RETRY_HPP_
