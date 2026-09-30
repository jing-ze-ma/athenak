//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file pgen_eos_utils.cpp
//! \brief the host helper of pgen_eos_utils.hpp that launches a kernel

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "pgen/pgen_eos_utils.hpp"

namespace pgen_eos {

Real HostGamma1FromP(const EOS_Data &eos, const Real d, const Real p) {
  if (!eos.IsGeneral()) return eos.gamma;
  DvceArray1D<Real> dout("pgen_eos_scalar", 1);
  par_for("pgen_eos_eval", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(int) {
    dout(0) = eos.Gamma1(d, eos.EnergyFromPressure(d, p));
  });
  auto hout = Kokkos::create_mirror_view(dout);
  Kokkos::deep_copy(hout, dout);
  return hout(0);
}

}  // namespace pgen_eos
