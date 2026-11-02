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

#include <string>

#include "parameter_input.hpp"

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

inline TracerLayout ReadTracerLayout(ParameterInput *pin) {
  TracerLayout t;
  if (pin->DoesParameterExist("mhd", "nscalars")) {
    t.nuser = pin->GetInteger("mhd", "nscalars");
  }
  t.xe = NamedTrue(pin, "mhd", "xe_tracer");
  t.fl = NamedTrue(pin, "mhd", "floor_tracer");
  return t;
}

}  // namespace mhd

#endif  // MHD_MHD_TRACERS_HPP_
