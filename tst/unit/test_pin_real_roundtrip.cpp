//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file test_pin_real_roundtrip.cpp
//  \brief standalone unit test: a Real recorded by GetOrAddReal()/SetReal() must come
//  back bit-identical after the parameter list has been dumped to text and re-parsed.
//  This is what happens across a restart, since the effective input is embedded in the
//  restart file.  See tst/unit/README for how to build and run it.

#include <Kokkos_Core.hpp>

#include <cmath>
#include <cstdio>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "parameter_input.hpp"

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);
  int nfail = 0;
  {
    std::vector<Real> vals = {
      static_cast<Real>(0.25*M_PI), static_cast<Real>(1.0/3.0),
      static_cast<Real>(0.1), static_cast<Real>(0.3),
      static_cast<Real>(6.02214076e23), static_cast<Real>(-2.7182818284590452),
      static_cast<Real>(1.0), static_cast<Real>(0.0),
      std::numeric_limits<Real>::min(), std::numeric_limits<Real>::denorm_min(),
      std::numeric_limits<Real>::max(), static_cast<Real>(1.0e-300)
    };
    // 1e-300 underflows to zero in single precision; drop it there
    ParameterInput pin;
    for (std::size_t n = 0; n < vals.size(); ++n) {
      std::stringstream nm;
      nm << "p" << n;
      Real got = pin.GetOrAddReal("problem", nm.str(), vals[n]);
      if (got != vals[n]) {
        std::printf("FAIL (in-memory) p%zu\n", n);
        nfail++;
      }
    }
    pin.SetReal("problem", "set0", static_cast<Real>(0.25*M_PI));

    // dump the whole parameter list to text, then re-parse it as a restart would
    std::stringstream dump;
    pin.ParameterDump(dump);
    ParameterInput pin2;
    pin2.LoadFromStream(dump);

    for (std::size_t n = 0; n < vals.size(); ++n) {
      std::stringstream nm;
      nm << "p" << n;
      Real back = pin2.GetReal("problem", nm.str());
      if (back != vals[n]) {
        nfail++;
        std::printf("FAIL p%zu: %.17g -> %.17g (text '%s')\n", n,
                    static_cast<double>(vals[n]), static_cast<double>(back),
                    ParameterInput::RealToString(vals[n]).c_str());
      }
    }
    if (pin2.GetReal("problem", "set0") != static_cast<Real>(0.25*M_PI)) {
      std::printf("FAIL set0 (SetReal)\n");
      nfail++;
    }
    // user-supplied text must be stored verbatim, not reformatted
    std::stringstream user;
    user << "<problem>\ncfl = 0.4\n";
    ParameterInput pin3;
    pin3.LoadFromStream(user);
    if (pin3.GetOrAddReal("problem", "cfl", static_cast<Real>(0.4)) !=
        static_cast<Real>(0.4)) {
      std::printf("FAIL verbatim user value\n");
      nfail++;
    }
    std::stringstream dump3;
    pin3.ParameterDump(dump3);
    if (dump3.str().find("0.4") == std::string::npos) {
      std::printf("FAIL user text '0.4' not preserved verbatim\n");
      nfail++;
    }
  }
  Kokkos::finalize();
  std::printf("%s (%d failures)\n", (nfail == 0) ? "PASS" : "FAIL", nfail);
  return (nfail == 0) ? 0 : 1;
}
