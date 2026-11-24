#ifndef RAD_VET_RAD_VET_ZALESAK_HPP_
#define RAD_VET_RAD_VET_ZALESAK_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_zalesak.hpp
//! \brief the cap warning of the Zalesak rounds of time2_face_fct and clean_wb_fct

#include <cstdio>

#include "globals.hpp"

namespace radvet {

//! \brief the Zalesak rounds stopped at their cap with R still falling: the limited E is
//! then not bounded (never silent)
inline void ZalesakCapWarn(const char *who, const int cycle, const int t2s) {
  if (global_variable::my_rank == 0) {
    std::printf("CLEAN WARNING cycle %d t2s %d: %s Zalesak rounds at the cap 4096, R not "
                "converged (E above the floor not guaranteed)\n", cycle, t2s, who);
  }
}

} // namespace radvet
#endif // RAD_VET_RAD_VET_ZALESAK_HPP_
