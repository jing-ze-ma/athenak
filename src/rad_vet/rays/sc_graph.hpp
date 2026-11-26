#ifndef RAD_VET_RAYS_SC_GRAPH_HPP_
#define RAD_VET_RAYS_SC_GRAPH_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_graph.hpp
//! \brief a launch sequence on one GPU stream captured as one CUDA / HIP graph and
//! replayed while its key (every input of the sequence's launches) is unchanged.
//! Start(key): a key seen with a graph replays it (2); a key seen once before is
//! captured (1: run the sequence, then End); a new key runs launch by launch (0) and is
//! remembered, so only a sequence that ran warm once (lazy allocations made) is
//! captured. The captured sequence must not synchronise the host with the stream.
//! Host(): per-key host values of the caller (e.g. counter increments of the captured
//! sequence, re-applied at each replay).

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "athena.hpp"

namespace radvet {

#if defined(KOKKOS_ENABLE_CUDA)
#define VET_GRAPH_ON 1
#define VGR(x) cuda##x
#define VGR_STREAM(ex) (ex).cuda_stream()
#elif defined(KOKKOS_ENABLE_HIP)
#define VET_GRAPH_ON 1
#define VGR(x) hip##x
#define VGR_STREAM(ex) (ex).hip_stream()
#else
#define VET_GRAPH_ON 0
#endif

class ScGraph {
 public:
  ScGraph() = default;
  ScGraph(const ScGraph &) = delete;
  ScGraph &operator=(const ScGraph &) = delete;
  ~ScGraph() {Clear();}
  // 2: replayed on ex's stream; 1: capture started (run the sequence, then End);
  // 0: run the sequence launch by launch
  int Start(const std::vector<int64_t> &key, const DevExeSpace &ex) {
#if VET_GRAPH_ON
    ++tick_;
    int q = -1;
    for (size_t s = 0; s < slot_.size(); ++s) {
      if (slot_[s].key == key) {q = static_cast<int>(s);}
    }
    if (q >= 0 && slot_[q].ex != nullptr) {
      slot_[q].use = tick_;
      last_ = q;
      Check(VGR(GraphLaunch)(slot_[q].ex, VGR_STREAM(ex)), "launch");
      nreplay_ += 1.0;
      return 2;
    }
    if (q >= 0) {
      slot_[q].use = tick_;
      Check(VGR(StreamBeginCapture)(VGR_STREAM(ex), VGR(StreamCaptureModeRelaxed)),
            "begin capture");
      cur_ = q;
      last_ = q;
      return 1;
    }
    // a new key (least recently used slot out)
    if (static_cast<int>(slot_.size()) < kSlots) {
      slot_.emplace_back();
      q = static_cast<int>(slot_.size()) - 1;
    } else {
      q = 0;
      for (int s = 1; s < kSlots; ++s) {if (slot_[s].use < slot_[q].use) {q = s;}}
      Drop(q);
    }
    slot_[q].key = key;
    slot_[q].use = tick_;
#else
    (void)key; (void)ex;
#endif
    return 0;
  }
  // end of a sequence started with mode 1: the graph is made and launched
  void End(const DevExeSpace &ex) {
#if VET_GRAPH_ON
    if (cur_ < 0) {return;}
    VGR(Graph_t) g = nullptr;
    Check(VGR(StreamEndCapture)(VGR_STREAM(ex), &g), "end capture");
    Check(VGR(GraphInstantiateWithFlags)(&slot_[cur_].ex, g, 0), "instantiate");
    (void)VGR(GraphDestroy)(g);
    Check(VGR(GraphLaunch)(slot_[cur_].ex, VGR_STREAM(ex)), "launch");
    ncapture_ += 1.0;
    cur_ = -1;
#else
    (void)ex;
#endif
  }
  // every graph out (a launch input changed that no key holds)
  void Clear() {
#if VET_GRAPH_ON
    for (size_t s = 0; s < slot_.size(); ++s) {Drop(static_cast<int>(s));}
    slot_.clear();
#endif
  }
  bool Capturing() const {
#if VET_GRAPH_ON
    return cur_ >= 0;
#else
    return false;
#endif
  }
  // host values of the key last replayed or captured
  std::vector<double> &Host() {
#if VET_GRAPH_ON
    return slot_[last_].host;
#else
    return host_;
#endif
  }
  double ncapture_ = 0.0, nreplay_ = 0.0;   // diagnostics

 private:
  static constexpr int kSlots = 4;
#if VET_GRAPH_ON
  struct Slot {
    std::vector<int64_t> key;
    VGR(GraphExec_t) ex = nullptr;
    int64_t use = 0;
    std::vector<double> host;
  };
  std::vector<Slot> slot_;
  int cur_ = -1, last_ = 0;
  int64_t tick_ = 0;
  void Drop(const int s) {
    if (slot_[s].ex != nullptr) {(void)VGR(GraphExecDestroy)(slot_[s].ex);}
    slot_[s].ex = nullptr;
    slot_[s].key.clear();
    slot_[s].host.clear();
  }
  static void Check(const VGR(Error_t) e, const char *what) {
    if (e != VGR(Success)) {
      std::cout << "### FATAL ERROR in sc_graph.hpp" << std::endl
                << "<rad_vet>/vet_gd graph: " << what << " failed: "
                << VGR(GetErrorString)(e) << std::endl;
      std::_Exit(EXIT_FAILURE);
    }
  }
#else
  std::vector<double> host_;
#endif
};

}  // namespace radvet

#endif  // RAD_VET_RAYS_SC_GRAPH_HPP_
