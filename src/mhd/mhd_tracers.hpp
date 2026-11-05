#ifndef MHD_MHD_TRACERS_HPP_
#define MHD_MHD_TRACERS_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_tracers.hpp
//! \brief the MHD passive-tracer slot layout (tracer-1101), computed in ONE place for the
//! MHD constructor, the RKG registers (Resistivity), the restart reader (pgen.cpp) and
//! main.cpp (the layout the restart file was written with).
//!
//! Slot order in u0/w0 (fixed):
//!   [0, nmhd)                        the MHD variables
//!   [nmhd, nmhd+nuser)               the user's <mhd>/nscalars
//!   nmhd+nuser                       <mhd>/xe_tracer: y = n_e/n_H (if on)
//!   last                             <mhd>/floor_tracer: floor-added mass (if on)
//! Both tracer keys are read only when named, so an input without them (and its
//! restarts) has nscalars = <mhd>/nscalars exactly as before.
//!
//! share-1104 P4: the SAME layout serves <hydro> (ReadTracerLayout(pin, "hydro")):
//! [0, nhydro) hydro variables, then <hydro>/nscalars, then <hydro>/floor_tracer (last).
//! Hydro has no x_e tracer (its consumer, the resistivity, is MHD only: refused there).
//! The floor-tracer drain and the x1-ghost rule are shared too (fluid:: below,
//! mhd_tracers.cpp); hydro's floor increments are added in the hydro ConsToPrim kernels.

#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"

class MeshBlockPack;

namespace mhd {

struct TracerLayout {
  int nuser = 0;        // <mhd>/nscalars
  bool xe = false;      // <mhd>/xe_tracer
  bool fl = false;      // <mhd>/floor_tracer
  int Total() const {return nuser + (xe ? 1 : 0) + (fl ? 1 : 0);}
  int NTracer() const {return (xe ? 1 : 0) + (fl ? 1 : 0);}
  // offsets from nmhd (-1 = off)
  int XeOff() const {return xe ? nuser : -1;}
  int FlOff() const {return fl ? nuser + (xe ? 1 : 0) : -1;}
  // compact code for globals (bit 0 = xe, bit 1 = floor)
  int Mask() const {return (xe ? 1 : 0) | (fl ? 2 : 0);}
};

// a boolean key that is named AND true; nothing is added to the input dump
inline bool NamedTrue(ParameterInput *pin, const std::string &blk,
                      const std::string &key) {
  return pin->DoesParameterExist(blk, key) && pin->GetBoolean(blk, key);
}

inline TracerLayout ReadTracerLayout(ParameterInput *pin,
                                     const std::string &blk = "mhd") {
  TracerLayout t;
  if (pin->DoesParameterExist(blk, "nscalars")) {
    t.nuser = pin->GetInteger(blk, "nscalars");
  }
  t.xe = (blk == "mhd") && NamedTrue(pin, blk, "xe_tracer");
  t.fl = NamedTrue(pin, blk, "floor_tracer");
  return t;
}

}  // namespace mhd

namespace fluid {
//! share-1104 P4: the floor-tracer drain kernel of MHD::FloorTracerDrain, shared with
//! Hydro::FloorTracerDrain.  bcc empty (hydro): the floor density is dfloor alone;
//! otherwise max(dfloor, B^2 va2i) as before.  phi empty unless etotgrav.  Adds the
//! active-cell drained mass/energy (volume integrated) to dm_sum/de_sum.
void FloorTracerDrainCore(MeshBlockPack *pp, const DvceArray5D<Real> &u0,
                          const DvceArray5D<Real> &w0, const DvceArray5D<Real> &bcc,
                          const DvceArray5D<Real> &wder, const DvceArray4D<Real> &phi,
                          const int isl, const int ns0, const int ns1, const Real dfl,
                          const Real va2i, const Real marg, const Real fmx,
                          Real &dm_sum, Real &de_sum);
//! the tracer slots [t0, t1] of the radial ghosts at physical x1 boundaries take the
//! edge cell's specific value (MHD::TracerX1Ghosts / Hydro::TracerX1Ghosts)
void TracerX1GhostsCore(MeshBlockPack *pp, const DvceArray5D<Real> &u0,
                        const DvceArray5D<Real> &w0, const int t0, const int t1,
                        const bool set_w);
}  // namespace fluid

#endif  // MHD_MHD_TRACERS_HPP_
