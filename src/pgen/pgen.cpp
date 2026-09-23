//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file pgen.cpp
//! \brief Implementation of constructors and functions in class ProblemGenerator.
//! Default constructor calls problem generator function, while  constructor for restarts
//! reads data from restart file, as well as re-initializing problem-specific data.

#include <iostream>
#include <string>
#include <utility>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

#include "athena.hpp"
#include "geodesic-grid/geodesic_grid.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "coordinates/adm.hpp"
#include "z4c/compact_object_tracker.hpp"
#include "z4c/z4c.hpp"
#include "radiation/radiation.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"
#include "srcterms/turb_driver.hpp"
#include "utils/two_stream_warm_rst.hpp"
#include "pgen.hpp"


//----------------------------------------------------------------------------------------
// default constructor, calls pgen function.

ProblemGenerator::ProblemGenerator(ParameterInput *pin, Mesh *pm) :
    hot_jupiter(false),
    user_bcs(false),
    user_srcs(false),
    user_hist(false),
    pmy_mesh_(pm) {
        
  hot_jupiter = pin->GetOrAddBoolean("problem","hot_jupiter",false);
  if (hot_jupiter) {
      hot_jupiter_param.Teq = pin->GetReal("problem","Teq");
      hot_jupiter_param.omega = pin->GetReal("problem","omega");
      hot_jupiter_param.grav = pin->GetReal("problem","grav");
      hot_jupiter_param.ap = pin->GetReal("problem","ap");
      hot_jupiter_param.Rgas = pin->GetReal("problem","Rgas");
      hot_jupiter_param.met = pin->GetReal("problem","met");
      hot_jupiter_param.bbot = pin->GetReal("problem","bbot");
      hot_jupiter_param.grav_point_mass =
          pin->GetOrAddBoolean("problem","grav_point_mass",false);
      hot_jupiter_param.stellar_tide =
          pin->GetOrAddBoolean("problem","stellar_tide",false);
      hot_jupiter_param.rot_potential =
          pin->GetOrAddBoolean("problem","rot_potential",false);
  }
  // check for user-defined boundary conditions
  for (int dir=0; dir<6; ++dir) {
    if (pm->mesh_bcs[dir] == BoundaryFlag::user) {
      user_bcs = true;
    }
  }

  user_srcs = pin->GetOrAddBoolean("problem","user_srcs",false);
  user_hist = pin->GetOrAddBoolean("problem","user_hist",false);

  // second argument false since this IS NOT a restart
  CallProblemGenerator(pin, false);

  // Check that user defined BCs were enrolled if needed
  if (user_bcs) {
    if (user_bcs_func == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "User BCs specified in <mesh> block, but not enrolled "
                << "by SetProblemData()." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  // Check that user defined srcterms were enrolled if needed
  if (user_srcs) {
    if (user_srcs_func == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "User SRCs specified in <problem> block, but not "
                << "enrolled by UserProblem()." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  // Check that user defined history outputs were enrolled if needed
  if (user_hist) {
    if (user_hist_func == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "User history output specified in <problem> block, but "
                << "not enrolled by UserProblem()." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
}

//----------------------------------------------------------------------------------------
// constructor for restarts
// When called, data needed to rebuild mesh has been read from restart file by
// Mesh::BuildTreeFromRestart() function. This constructor reads from the restart file and
// initializes all the dependent variables (u0,b0,etc) stored in each Physics class. It
// also calls ProblemGenerator::SetProblemData() function to set any user-defined BCs,
// and any data necessary for restart runs to continue correctly.

ProblemGenerator::ProblemGenerator(ParameterInput *pin, Mesh *pm, IOWrapper resfile,
                                   bool single_file_per_rank) :
    hot_jupiter(false),
    user_bcs(false),
    user_srcs(false),
    user_hist(false),
    pmy_mesh_(pm) {
        
  hot_jupiter = pin->GetOrAddBoolean("problem","hot_jupiter",false);
  if (hot_jupiter) {
    hot_jupiter_param.Teq = pin->GetReal("problem","Teq");
    hot_jupiter_param.omega = pin->GetReal("problem","omega");
    hot_jupiter_param.grav = pin->GetReal("problem","grav");
    hot_jupiter_param.ap = pin->GetReal("problem","ap");
    hot_jupiter_param.Rgas = pin->GetReal("problem","Rgas");
    hot_jupiter_param.met = pin->GetReal("problem","met");
    hot_jupiter_param.bbot = pin->GetReal("problem","bbot");
    hot_jupiter_param.grav_point_mass =
        pin->GetOrAddBoolean("problem","grav_point_mass",false);
    hot_jupiter_param.stellar_tide =
        pin->GetOrAddBoolean("problem","stellar_tide",false);
    hot_jupiter_param.rot_potential =
        pin->GetOrAddBoolean("problem","rot_potential",false);
  }
        
  // check for user-defined boundary conditions
  for (int dir=0; dir<6; ++dir) {
    if (pm->mesh_bcs[dir] == BoundaryFlag::user) {
      user_bcs = true;
    }
  }
  user_srcs = pin->GetOrAddBoolean("problem","user_srcs",false);
  user_hist = pin->GetOrAddBoolean("problem","user_hist",false);

  // get spatial dimensions of arrays, including ghost zones
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int nout1 = indcs.nx1 + 2*(indcs.ng);
  int nout2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int nout3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  int nmb = pm->pmb_pack->nmb_thispack;
  // calculate total number of CC variables
  hydro::Hydro* phydro = pm->pmb_pack->phydro;
  mhd::MHD* pmhd = pm->pmb_pack->pmhd;
  adm::ADM* padm = pm->pmb_pack->padm;
  z4c::Z4c* pz4c = pm->pmb_pack->pz4c;
  radiation::Radiation* prad=pm->pmb_pack->prad;
  TurbulenceDriver* pturb=pm->pmb_pack->pturb;
  // <rad_m1>: the four evolved moments, read back in the SAME order in which
  // restart.cpp writes them (after the <radiation> intensities, before the turbulence
  // forcing).  Guarded by the pointer, so a file from a run without <rad_m1> is
  // unchanged and still readable.
  radm1::RadiationM1* pradm1 = pm->pmb_pack->pradm1;
  int nrad = 0, nm1 = 0, nhydro = 0, nmhd = 0, nforce = 3, nadm = 0, nz4c = 0;
  if (pradm1 != nullptr) {
    nm1 = radm1::M1_NVAR;
  }
  if (phydro != nullptr) {
    nhydro = phydro->nhydro + phydro->nscalars;
  }
  if (pmhd != nullptr) {
    nmhd = pmhd->nmhd + pmhd->nscalars;
  }
  if (prad != nullptr) {
    nrad = prad->prgeo->nangles;
  }
  if (pz4c != nullptr) {
    nz4c = pz4c->nz4c;
  } else if (padm != nullptr) {
    nadm = padm->nadm;
  }

  // root process reads z4c last_output_time and tracker data
  if (pz4c != nullptr) {
    Real last_output_time;
    if (global_variable::my_rank == 0 || single_file_per_rank) {
      if (resfile.Read_Reals(&last_output_time, 1,single_file_per_rank) != 1) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "z4c::last_output_time data size read from restart "
                  << "file is incorrect, restart file is broken." << std::endl;
        exit(EXIT_FAILURE);
      }
    }
#if MPI_PARALLEL_ENABLED
    if (!single_file_per_rank) {
      MPI_Bcast(&last_output_time, sizeof(Real), MPI_CHAR, 0, MPI_COMM_WORLD);
    }
#endif
    pz4c->last_output_time = last_output_time;

    for (auto &pt : pz4c->ptracker) {
      Real pos[3];
      if (global_variable::my_rank == 0 || single_file_per_rank) {
        if (resfile.Read_Reals(&pos[0], 3, single_file_per_rank) != 3) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "compact object tracker data size read from restart "
                    << "file is incorrect, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
      }
#if MPI_PARALLEL_ENABLED
      if (!single_file_per_rank) {
        MPI_Bcast(&pos[0], 3*sizeof(Real), MPI_CHAR, 0, MPI_COMM_WORLD);
      }
#endif
      pt->SetPos(&pos[0]);
    }
  }

  if (pturb != nullptr) {
    // root process reads size the random seed
    char *rng_data = new char[sizeof(RNG_State)];
    // the master process reads the variables data
    if (global_variable::my_rank == 0 || single_file_per_rank) {
      if (resfile.Read_bytes(rng_data, 1, sizeof(RNG_State), single_file_per_rank)
          != sizeof(RNG_State)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "RNG data size read from restart file is incorrect, "
                  << "restart file is broken." << std::endl;
        exit(EXIT_FAILURE);
      }
    }
#if MPI_PARALLEL_ENABLED
    if (!single_file_per_rank) {
      // then broadcast the RNG information
      MPI_Bcast(rng_data, sizeof(RNG_State), MPI_CHAR, 0, MPI_COMM_WORLD);
    }
#endif
    std::memcpy(&(pturb->rstate), &(rng_data[0]), sizeof(RNG_State));
  }

  // --- THE OPTIONAL PGEN STATE BLOCK (see kPgenRstMagic in pgen.hpp).  It is marked, so
  // the eight bytes read here decide: the marker means a block follows and the variable
  // data size comes after it, anything else IS the variable data size of the old layout
  // and has already been read.  That is why the read is unconditional and needs no switch
  // from the input file -- a file written without the block, which is every file every
  // other problem generator writes and every file written before this existed, takes the
  // second branch and sees exactly the bytes it always did.
  IOWrapperSizeT variablesize = sizeof(IOWrapperSizeT);
  char *variabledata = new char[variablesize];
  static_assert(sizeof(kPgenRstMagic) == sizeof(IOWrapperSizeT),
                "the pgen state marker must be exactly as long as the data size it "
                "stands in front of, or the peek below cannot fall back");
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    if (resfile.Read_bytes(variabledata, 1, variablesize, single_file_per_rank)
        != variablesize) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Variable data size read from restart file is incorrect, "
                << "restart file is broken." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
#if MPI_PARALLEL_ENABLED
  // then broadcast the datasize information
  if (!single_file_per_rank) {
    MPI_Bcast(variabledata, variablesize, MPI_CHAR, 0, MPI_COMM_WORLD);
  }
#endif
  if (std::memcmp(variabledata, &(kPgenRstMagic[0]), sizeof(kPgenRstMagic)) == 0) {
    // the marker: read the length, then the payload, then the real data size
    IOWrapperSizeT nb = 0;
    if (global_variable::my_rank == 0 || single_file_per_rank) {
      if (resfile.Read_bytes(&nb, 1, sizeof(IOWrapperSizeT), single_file_per_rank)
          != sizeof(IOWrapperSizeT)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "the problem generator state block of this restart "
                  << "file has no length, restart file is broken." << std::endl;
        exit(EXIT_FAILURE);
      }
    }
#if MPI_PARALLEL_ENABLED
    if (!single_file_per_rank) {
      MPI_Bcast(&nb, sizeof(IOWrapperSizeT), MPI_CHAR, 0, MPI_COMM_WORLD);
    }
#endif
    pgen_rststate.resize(nb);
    if (nb > 0) {
      if (global_variable::my_rank == 0 || single_file_per_rank) {
        if (resfile.Read_bytes(pgen_rststate.data(), 1, nb, single_file_per_rank) != nb) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "the problem generator state block of this restart "
                    << "file is short, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
      }
#if MPI_PARALLEL_ENABLED
      if (!single_file_per_rank) {
        MPI_Bcast(pgen_rststate.data(), nb, MPI_CHAR, 0, MPI_COMM_WORLD);
      }
#endif
    }
    if (global_variable::my_rank == 0 || single_file_per_rank) {
      if (resfile.Read_bytes(variabledata, 1, variablesize, single_file_per_rank)
          != variablesize) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Variable data size read from restart file is "
                  << "incorrect, restart file is broken." << std::endl;
        exit(EXIT_FAILURE);
      }
    }
#if MPI_PARALLEL_ENABLED
    if (!single_file_per_rank) {
      MPI_Bcast(variabledata, variablesize, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
#endif
  }

  // --- THE MODE-3 WARM-START HEADER (kRtWarmRstMagic, utils/two_stream_warm_rst.hpp).
  // The same eight-byte peek once more, immediately behind the pgen state block: the
  // marker means the little header follows and the variable data size comes after it,
  // anything else IS the variable data size and has already been read.  The header says
  // how many per-cell history levels the tail carries -- which the tail's LENGTH cannot,
  // since one extra slab is also what an mhd wtemp looks like and two are what the hydro
  // derived cache looks like -- and carries the previous call's bdt with them.
  int nwarm_file = 0;
  Real warm_bdt_prev = 0.0;
  if (std::memcmp(variabledata, &(two_stream_rt::kRtWarmRstMagic[0]),
                  sizeof(two_stream_rt::kRtWarmRstMagic)) == 0) {
    char warm_hdr[2*sizeof(std::int32_t) + sizeof(Real)];
    IOWrapperSizeT nb = 0;
    if (global_variable::my_rank == 0 || single_file_per_rank) {
      if (resfile.Read_bytes(&nb, 1, sizeof(IOWrapperSizeT), single_file_per_rank)
          != sizeof(IOWrapperSizeT)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "the warm-start header of this restart file has no "
                  << "length, restart file is broken." << std::endl;
        exit(EXIT_FAILURE);
      }
    }
#if MPI_PARALLEL_ENABLED
    if (!single_file_per_rank) {
      MPI_Bcast(&nb, sizeof(IOWrapperSizeT), MPI_CHAR, 0, MPI_COMM_WORLD);
    }
#endif
    if (nb != sizeof(warm_hdr)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "the warm-start header of this restart file is "
                << nb << " bytes, not " << sizeof(warm_hdr) << ", restart file is "
                << "broken." << std::endl;
      exit(EXIT_FAILURE);
    }
    if (global_variable::my_rank == 0 || single_file_per_rank) {
      if (resfile.Read_bytes(&(warm_hdr[0]), 1, nb, single_file_per_rank) != nb) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "the warm-start header of this restart file is short, "
                  << "restart file is broken." << std::endl;
        exit(EXIT_FAILURE);
      }
      if (resfile.Read_bytes(variabledata, 1, variablesize, single_file_per_rank)
          != variablesize) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Variable data size read from restart file is "
                  << "incorrect, restart file is broken." << std::endl;
        exit(EXIT_FAILURE);
      }
    }
#if MPI_PARALLEL_ENABLED
    if (!single_file_per_rank) {
      MPI_Bcast(&(warm_hdr[0]), nb, MPI_CHAR, 0, MPI_COMM_WORLD);
      MPI_Bcast(variabledata, variablesize, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
#endif
    std::int32_t hdr[2];
    std::memcpy(&(hdr[0]), &(warm_hdr[0]), sizeof(hdr));
    std::memcpy(&warm_bdt_prev, &(warm_hdr[0]) + sizeof(hdr), sizeof(warm_bdt_prev));
    nwarm_file = static_cast<int>(hdr[0]);
    if (nwarm_file < 0 || nwarm_file > 2) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "the warm-start header of this restart file claims "
                << nwarm_file << " history levels, restart file is broken." << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // --- THE <rad_m1> implicit_predictor HEADER (radm1::kM1PredRstMagic): the same peek
  // once more, behind the warm-start header.  It says how many ipred channels the tail
  // carries behind the warm-start slabs and the dt of the step they were stored on.
  int npred_file = 0;
  Real pred_dt_file = 0.0;
  if (std::memcmp(variabledata, &(radm1::kM1PredRstMagic[0]),
                  sizeof(radm1::kM1PredRstMagic)) == 0) {
    char pred_hdr[2*sizeof(std::int32_t) + sizeof(Real)];
    IOWrapperSizeT nb = 0;
    bool ok = true;
    if (global_variable::my_rank == 0 || single_file_per_rank) {
      ok = (resfile.Read_bytes(&nb, 1, sizeof(IOWrapperSizeT), single_file_per_rank)
            == sizeof(IOWrapperSizeT)) && (nb == sizeof(pred_hdr));
      ok = ok && (resfile.Read_bytes(&(pred_hdr[0]), 1, nb, single_file_per_rank) == nb);
      ok = ok && (resfile.Read_bytes(variabledata, 1, variablesize, single_file_per_rank)
                  == variablesize);
    }
#if MPI_PARALLEL_ENABLED
    if (!single_file_per_rank) {
      MPI_Bcast(&ok, sizeof(bool), MPI_CHAR, 0, MPI_COMM_WORLD);
      MPI_Bcast(&(pred_hdr[0]), sizeof(pred_hdr), MPI_CHAR, 0, MPI_COMM_WORLD);
      MPI_Bcast(variabledata, variablesize, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
#endif
    std::int32_t hdr[2] = {0, 0};
    if (ok) {
      std::memcpy(&(hdr[0]), &(pred_hdr[0]), sizeof(hdr));
      std::memcpy(&pred_dt_file, &(pred_hdr[0]) + sizeof(hdr), sizeof(pred_dt_file));
    }
    if (!ok || hdr[0] != 1 || hdr[1] < 1 || hdr[1] > 3) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "the <rad_m1> predictor header of this restart file is "
                << "broken." << std::endl;
      exit(EXIT_FAILURE);
    }
    npred_file = static_cast<int>(hdr[1]);
  }

  // --- THE <rad_m1> time_scheme = hesdirk2 HEADER (radm1::kM1Time2RstMagic), behind the
  // predictor one: nch slabs of the stored FSAL slope follow the predictor slabs.
  int nt2_file = 0;
  Real t2_hv[3] = {0.0, 0.0, 0.0};
  if (std::memcmp(variabledata, &(radm1::kM1Time2RstMagic[0]),
                  sizeof(radm1::kM1Time2RstMagic)) == 0) {
    char t2_hdr[2*sizeof(std::int32_t) + 3*sizeof(Real)];
    IOWrapperSizeT nb = 0;
    bool ok = true;
    if (global_variable::my_rank == 0 || single_file_per_rank) {
      ok = (resfile.Read_bytes(&nb, 1, sizeof(IOWrapperSizeT), single_file_per_rank)
            == sizeof(IOWrapperSizeT)) && (nb == sizeof(t2_hdr));
      ok = ok && (resfile.Read_bytes(&(t2_hdr[0]), 1, nb, single_file_per_rank) == nb);
      ok = ok && (resfile.Read_bytes(variabledata, 1, variablesize, single_file_per_rank)
                  == variablesize);
    }
#if MPI_PARALLEL_ENABLED
    if (!single_file_per_rank) {
      MPI_Bcast(&ok, sizeof(bool), MPI_CHAR, 0, MPI_COMM_WORLD);
      MPI_Bcast(&(t2_hdr[0]), sizeof(t2_hdr), MPI_CHAR, 0, MPI_COMM_WORLD);
      MPI_Bcast(variabledata, variablesize, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
#endif
    std::int32_t hdr[2] = {0, 0};
    if (ok) {
      std::memcpy(&(hdr[0]), &(t2_hdr[0]), sizeof(hdr));
      std::memcpy(&(t2_hv[0]), &(t2_hdr[0]) + sizeof(hdr), sizeof(t2_hv));
    }
    if (!ok || hdr[0] != 1 || hdr[1] < 1 || hdr[1] > 64) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "the <rad_m1> hesdirk2 header of this restart file is "
                << "broken." << std::endl;
      exit(EXIT_FAILURE);
    }
    nt2_file = static_cast<int>(hdr[1]);
  }

  IOWrapperSizeT data_size;
  std::memcpy(&data_size, &(variabledata[0]), sizeof(IOWrapperSizeT));

  // calculate total number of CC variables
  IOWrapperSizeT headeroffset;
  // master process gets file offset
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    headeroffset = resfile.GetPosition(single_file_per_rank);
  }
#if MPI_PARALLEL_ENABLED
  // then broadcasts it
  if (!single_file_per_rank) {
    MPI_Bcast(&headeroffset, sizeof(IOWrapperSizeT), MPI_CHAR, 0, MPI_COMM_WORLD);
  }
#endif

  IOWrapperSizeT data_size_ = 0;
  if (phydro != nullptr) {
    data_size_ += nout1*nout2*nout3*nhydro*sizeof(Real); // hydro u0
  }
  if (pmhd != nullptr) {
    data_size_ += nout1*nout2*nout3*nmhd*sizeof(Real);   // mhd u0
    data_size_ += (nout1+1)*nout2*nout3*sizeof(Real);    // mhd b0.x1f
    data_size_ += nout1*(nout2+1)*nout3*sizeof(Real);    // mhd b0.x2f
    data_size_ += nout1*nout2*(nout3+1)*sizeof(Real);    // mhd b0.x3f
  }
  if (prad != nullptr) {
    data_size_ += nout1*nout2*nout3*nrad*sizeof(Real);   // rad i0
  }
  if (pradm1 != nullptr) {
    data_size_ += nout1*nout2*nout3*nm1*sizeof(Real);    // rad_m1 u0
    if (pradm1->transport >= radm1::M1_TRANSPORT_IMPLICIT_X1) {
      data_size_ += (nout1+1)*nout2*nout3*sizeof(Real);  // rad_m1 f0x1 (milestone 3a)
    }
  }
  // milestone 3b phase B: the x2/x3 face-normal fluxes of transport = implicit.  A file
  // written before they existed (or by an implicit_x1 run) simply does not have them:
  // that length is accepted below, with one warning, and the arrays stay zero.
  IOWrapperSizeT m1f23_size = 0;
  if (pradm1 != nullptr && pradm1->trans_on) {
    m1f23_size += nout1*(nout2+1)*nout3*sizeof(Real);
    if (pradm1->trans_x3) {m1f23_size += nout1*nout2*(nout3+1)*sizeof(Real);}
  }
  bool m1_have_f23 = (m1f23_size > 0);
  data_size_ += m1f23_size;
  if (pturb != nullptr) {
    data_size_ += nout1*nout2*nout3*nforce*sizeof(Real); // forcing
  }
  if (pz4c != nullptr) {
    data_size_ += nout1*nout2*nout3*nz4c*sizeof(Real);   // z4c u0
  } else if (padm != nullptr) {
    data_size_ += nout1*nout2*nout3*nadm*sizeof(Real);   // adm u_adm
  }

  // The general-EOS temperature cache (Hydro/MHD::wtemp) is appended to the END of each
  // MeshBlock record, because it is the warm start of the c2p temperature root find and
  // a cold start lands on a different temperature in the last digits -- a restart is then
  // not a bitwise continuation.  A file written before this was added (or by a run whose
  // EOS is ideal) simply does not have the tail: accept that size too, and leave wtemp
  // as it is, which costs one cold start on the first conversion and nothing after.
  //
  // THE DERIVED CACHE (Hydro::wder and MHD::wder, p and Gamma_1) follows it, for
  // the reason spelled out on outarray_wdp in outputs.hpp: restoring wtemp alone is not
  // enough, because the first Fluxes call of the restarted cycle reconstructs p and
  // Gamma_1 that no ConsToPrim of the restarted run has produced yet, and re-deriving
  // them from the restored temperature moves them by an ULP (the log10/Pow10 round trip
  // of T is not the identity).  The hydro derived cache comes first, the MHD one behind
  // it.  Four tail lengths are therefore accepted -- both caches; the hydro derived
  // cache but not the MHD one (a file written after the hydro feature and before the MHD
  // one, which only a hydro+MHD ion-neutral run can see); wtemp only (a file written
  // before either derived cache existed, which is what every general-EOS MHD file up to
  // now is); and neither -- and only the first makes the restart a bitwise continuation.
  bool wt_hyd = (phydro != nullptr) && phydro->peos->eos_data.IsGeneral();
  bool wt_mhd = (pmhd != nullptr) && pmhd->peos->eos_data.IsGeneral();
  bool wd_hyd = wt_hyd;
  bool wd_mhd = wt_mhd;
  IOWrapperSizeT wt_size = 0;
  if (wt_hyd) { wt_size += nout1*nout2*nout3*sizeof(Real); }
  if (wt_mhd) { wt_size += nout1*nout2*nout3*sizeof(Real); }
  IOWrapperSizeT wd_size = 0;
  if (wd_hyd) { wd_size += 2*nout1*nout2*nout3*sizeof(Real); }
  if (wd_mhd) { wd_size += 2*nout1*nout2*nout3*sizeof(Real); }
  // and behind both of them, the mode-3 warm-start history: nwarm_file slabs, a number
  // the marked header above gave us rather than something inferred from the length
  // (and behind those the npred_file <rad_m1> predictor slabs, also header-declared)
  IOWrapperSizeT wm_size = (nwarm_file + npred_file + nt2_file)*nout1*nout2*nout3
                           *sizeof(Real);
  if ((data_size_ + wt_size + wd_size + wm_size) == data_size) {
    data_size_ += wt_size + wd_size + wm_size;
  } else if (wd_hyd && wd_mhd &&
             (data_size_ + wt_size + wd_size/2 + wm_size) == data_size) {
    // a run with BOTH a general-EOS hydro and a general-EOS MHD module (ion-neutral)
    // reading a file written while only the hydro derived cache existed: the hydro half
    // of the tail is there, the MHD half is not
    data_size_ += wt_size + wd_size/2 + wm_size;
    wd_mhd = false;
    if (global_variable::my_rank == 0) {
      std::cout << "### WARNING: restart file has no general-EOS MHD derived "
                << "(p, Gamma_1) cache (written before it was added); the first flux "
                << "calculation re-derives it and this restart is not bitwise."
                << std::endl;
    }
  } else if (wd_size > 0 && (data_size_ + wt_size + wm_size) == data_size) {
    data_size_ += wt_size + wm_size;
    wd_hyd = false;
    wd_mhd = false;
    if (global_variable::my_rank == 0) {
      std::cout << "### WARNING: restart file has no general-EOS derived (p, Gamma_1) "
                << "cache (written before it was added); the first flux calculation "
                << "re-derives it and this restart is not bitwise." << std::endl;
    }
  } else {
    if (wt_size > 0 && (data_size_ + wm_size) == data_size) {
      if (global_variable::my_rank == 0) {
        std::cout << "### WARNING: restart file has no general-EOS temperature cache "
                  << "(written before it was added); the first conversion to primitives "
                  << "cold starts and this restart is not bitwise." << std::endl;
      }
    }
    wt_hyd = false;
    wt_mhd = false;
    wd_hyd = false;
    wd_mhd = false;
    data_size_ += wm_size;
  }

  // the file predates the transverse face fluxes (or was written by an implicit_x1 run)
  if (data_size_ != data_size && m1f23_size > 0 &&
      (data_size_ - m1f23_size) == data_size) {
    data_size_ -= m1f23_size;
    m1_have_f23 = false;
    if (global_variable::my_rank == 0) {
      std::cout << "### WARNING: restart file carries no <rad_m1> x2/x3 face fluxes "
                << "(written before transport = implicit existed); they are "
                << "zero-initialised and the first step rebuilds them." << std::endl;
    }
  }

  if (data_size_ != data_size) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "CC data size read from restart file not equal to size "
              << "of Hydro, MHD, Rad, and/or Z4c arrays, restart file is broken."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  // read CC data into host array
  int mygids = pm->gids_eachrank[global_variable::my_rank];
  IOWrapperSizeT offset_myrank = headeroffset;
  if (!single_file_per_rank) {
    offset_myrank += data_size_ * pm->gids_eachrank[global_variable::my_rank];
  }
  IOWrapperSizeT myoffset = offset_myrank;

  HostArray5D<Real> ccin("rst-cc-in", 1, 1, 1, 1, 1);
  HostFaceFld4D<Real> fcin("rst-fc-in", 1, 1, 1, 1);

  // calculate max/min number of MeshBlocks across all ranks
  int noutmbs_max = pm->nmb_eachrank[0];
  int noutmbs_min = pm->nmb_eachrank[0];
  for (int i=0; i<(global_variable::nranks); ++i) {
    noutmbs_max = std::max(noutmbs_max,pm->nmb_eachrank[i]);
    noutmbs_min = std::min(noutmbs_min,pm->nmb_eachrank[i]);
  }

  if (phydro != nullptr) {
    Kokkos::realloc(ccin, nmb, nhydro, nout3, nout2, nout1);
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to read, so read collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Read_Reals_at_all(mbptr.data(), mbcnt, myoffset, single_file_per_rank)
            != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC hydro data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Read_Reals_at(mbptr.data(), mbcnt, myoffset, single_file_per_rank)
            != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC hydro data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    Kokkos::deep_copy(Kokkos::subview(phydro->u0, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
    offset_myrank += nout1*nout2*nout3*nhydro*sizeof(Real); // hydro u0
    myoffset = offset_myrank;
  }

  if (pmhd != nullptr) {
    Kokkos::realloc(ccin, nmb, nmhd, nout3, nout2, nout1);
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to read, so read collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                   Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Read_Reals_at_all(mbptr.data(), mbcnt, myoffset, single_file_per_rank)
            != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC mhd data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Read_Reals_at(mbptr.data(), mbcnt, myoffset, single_file_per_rank)
            != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC mhd data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    Kokkos::deep_copy(Kokkos::subview(pmhd->u0, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
    offset_myrank += nout1*nout2*nout3*nmhd*sizeof(Real);   // mhd u0
    myoffset = offset_myrank;

    Kokkos::realloc(fcin.x1f, nmb, nout3, nout2, nout1+1);
    Kokkos::realloc(fcin.x2f, nmb, nout3, nout2+1, nout1);
    Kokkos::realloc(fcin.x3f, nmb, nout3+1, nout2, nout1);
    // read FC data into host array, again one MeshBlock at a time
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to x1-face field
        auto x1fptr = Kokkos::subview(fcin.x1f, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        int fldcnt = x1fptr.size();

        if (resfile.Read_Reals_at_all(x1fptr.data(), fldcnt, myoffset,
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Input b0.x1f field not read correctly from rst file, "
                << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        // get ptr to x2-face field
        auto x2fptr = Kokkos::subview(fcin.x2f, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        fldcnt = x2fptr.size();

        if (resfile.Read_Reals_at_all(x2fptr.data(), fldcnt, myoffset,
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Input b0.x2f field not read correctly from rst file, "
                << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        // get ptr to x3-face field
        auto x3fptr = Kokkos::subview(fcin.x3f, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        fldcnt = x3fptr.size();

        if (resfile.Read_Reals_at_all(x3fptr.data(), fldcnt, myoffset,
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Input b0.x3f field not read correctly from rst file, "
                << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        myoffset += data_size-(x1fptr.size()+x2fptr.size()+x3fptr.size())*sizeof(Real);
      } else if (m < pm->nmb_thisrank) {
        // get ptr to x1-face field
        auto x1fptr = Kokkos::subview(fcin.x1f, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        int fldcnt = x1fptr.size();

        if (resfile.Read_Reals_at(x1fptr.data(), fldcnt, myoffset,
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Input b0.x1f field not read correctly from rst file, "
                << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        // get ptr to x2-face field
        auto x2fptr = Kokkos::subview(fcin.x2f, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        fldcnt = x2fptr.size();

        if (resfile.Read_Reals_at(x2fptr.data(), fldcnt, myoffset,
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Input b0.x2f field not read correctly from rst file, "
                << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        // get ptr to x3-face field
        auto x3fptr = Kokkos::subview(fcin.x3f, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        fldcnt = x3fptr.size();

        if (resfile.Read_Reals_at(x3fptr.data(), fldcnt, myoffset,
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Input b0.x3f field not read correctly from rst file, "
                << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        myoffset += data_size-(x1fptr.size()+x2fptr.size()+x3fptr.size())*sizeof(Real);
      }
    }
    Kokkos::deep_copy(Kokkos::subview(pmhd->b0.x1f, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL), fcin.x1f);
    Kokkos::deep_copy(Kokkos::subview(pmhd->b0.x2f, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL), fcin.x2f);
    Kokkos::deep_copy(Kokkos::subview(pmhd->b0.x3f, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL), fcin.x3f);
    offset_myrank += (nout1+1)*nout2*nout3*sizeof(Real);    // mhd b0.x1f
    offset_myrank += nout1*(nout2+1)*nout3*sizeof(Real);    // mhd b0.x2f
    offset_myrank += nout1*nout2*(nout3+1)*sizeof(Real);    // mhd b0.x3f
    myoffset = offset_myrank;
  }

  if (prad != nullptr) {
    Kokkos::realloc(ccin, nmb, nrad, nout3, nout2, nout1);
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to read, so read collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Read_Reals_at_all(mbptr.data(), mbcnt, myoffset,
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC rad data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Read_Reals_at(mbptr.data(), mbcnt, myoffset,
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC rad data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    Kokkos::deep_copy(Kokkos::subview(prad->i0, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
    offset_myrank += nout1*nout2*nout3*nrad*sizeof(Real);   // radiation i0
    myoffset = offset_myrank;
  }

  if (pradm1 != nullptr) {
    Kokkos::realloc(ccin, nmb, nm1, nout3, nout2, nout1);
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to read, so read collectively
      if (m < noutmbs_min) {
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Read_Reals_at_all(mbptr.data(), mbcnt, myoffset,
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC rad_m1 data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished reading, so use non-collective read
      } else if (m < pm->nmb_thisrank) {
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Read_Reals_at(mbptr.data(), mbcnt, myoffset,
                                  single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC rad_m1 data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    Kokkos::deep_copy(Kokkos::subview(pradm1->u0, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
    offset_myrank += nout1*nout2*nout3*nm1*sizeof(Real);    // rad_m1 u0
    myoffset = offset_myrank;

    // milestone 3a / 3b phase B: the persistent face fluxes, read back in the same order
    // the writer (src/outputs/restart.cpp) put them down: f0x1, then f0x2, then f0x3.
    auto rdface = [&](HostArray4D<Real> &arr) {
      for (int m=0;  m<noutmbs_max; ++m) {
        if (m < noutmbs_min) {
          auto mbptr = Kokkos::subview(arr, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
          int mbcnt = mbptr.size();
          if (resfile.Read_Reals_at_all(mbptr.data(), mbcnt, myoffset,
                                        single_file_per_rank) != mbcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "rad_m1 face data not read correctly from rst "
                      << "file, restart file is broken." << std::endl;
            exit(EXIT_FAILURE);
          }
          myoffset += data_size;
        } else if (m < pm->nmb_thisrank) {
          auto mbptr = Kokkos::subview(arr, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
          int mbcnt = mbptr.size();
          if (resfile.Read_Reals_at(mbptr.data(), mbcnt, myoffset,
                                    single_file_per_rank) != mbcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "rad_m1 face data not read correctly from rst "
                      << "file, restart file is broken." << std::endl;
            exit(EXIT_FAILURE);
          }
          myoffset += data_size;
        }
      }
    };
    if (pradm1->transport >= radm1::M1_TRANSPORT_IMPLICIT_X1) {
      HostArray4D<Real> fcin("m1fcin", nmb, nout3, nout2, nout1+1);
      for (int m=0;  m<noutmbs_max; ++m) {
        if (m < noutmbs_min) {
          auto mbptr = Kokkos::subview(fcin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
          int mbcnt = mbptr.size();
          if (resfile.Read_Reals_at_all(mbptr.data(), mbcnt, myoffset,
                                        single_file_per_rank) != mbcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "rad_m1 face data not read correctly from rst "
                      << "file, restart file is broken." << std::endl;
            exit(EXIT_FAILURE);
          }
          myoffset += data_size;
        } else if (m < pm->nmb_thisrank) {
          auto mbptr = Kokkos::subview(fcin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
          int mbcnt = mbptr.size();
          if (resfile.Read_Reals_at(mbptr.data(), mbcnt, myoffset,
                                    single_file_per_rank) != mbcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "rad_m1 face data not read correctly from rst "
                      << "file, restart file is broken." << std::endl;
            exit(EXIT_FAILURE);
          }
          myoffset += data_size;
        }
      }
      Kokkos::deep_copy(Kokkos::subview(pradm1->f0x1, std::make_pair(0,nmb), Kokkos::ALL,
                        Kokkos::ALL, Kokkos::ALL), fcin);
      offset_myrank += (nout1+1)*nout2*nout3*sizeof(Real);   // rad_m1 f0x1
      myoffset = offset_myrank;
    }
    if (m1_have_f23) {
      HostArray4D<Real> f2in("m1f2in", nmb, nout3, nout2+1, nout1);
      rdface(f2in);
      Kokkos::deep_copy(Kokkos::subview(pradm1->f0x2, std::make_pair(0,nmb), Kokkos::ALL,
                        Kokkos::ALL, Kokkos::ALL), f2in);
      offset_myrank += nout1*(nout2+1)*nout3*sizeof(Real);   // rad_m1 f0x2
      myoffset = offset_myrank;
      if (pradm1->trans_x3) {
        HostArray4D<Real> f3in("m1f3in", nmb, nout3+1, nout2, nout1);
        rdface(f3in);
        Kokkos::deep_copy(Kokkos::subview(pradm1->f0x3, std::make_pair(0,nmb),
                          Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), f3in);
        offset_myrank += nout1*nout2*(nout3+1)*sizeof(Real);   // rad_m1 f0x3
        myoffset = offset_myrank;
      }
    }
  }

  if (pturb != nullptr) {
    Kokkos::realloc(ccin, nmb, nforce, nout3, nout2, nout1);
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to read, so read collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Read_Reals_at_all(mbptr.data(), mbcnt, myoffset,
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC turb data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Read_Reals_at(mbptr.data(), mbcnt, myoffset,
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC turb data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    Kokkos::deep_copy(Kokkos::subview(pturb->force, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
    offset_myrank += nout1*nout2*nout3*nforce*sizeof(Real); // forcing
    myoffset = offset_myrank;
  }

  if (pz4c != nullptr) {
    Kokkos::realloc(ccin, nmb, nz4c, nout3, nout2, nout1);
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to read, so read collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Read_Reals_at_all(mbptr.data(), mbcnt, myoffset,
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC z4c data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Read_Reals_at(mbptr.data(), mbcnt, myoffset,
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC z4c data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    Kokkos::deep_copy(Kokkos::subview(pz4c->u0, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
    offset_myrank += nout1*nout2*nout3*nz4c*sizeof(Real);   // z4c u0
    myoffset = offset_myrank;

    // We also need to reinitialize the ADM data.
    pz4c->Z4cToADM(pmy_mesh_->pmb_pack);
  } else if (padm != nullptr) {
    Kokkos::realloc(ccin, nmb, nadm, nout3, nout2, nout1);
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to read, so read collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Read_Reals_at_all(mbptr.data(), mbcnt, myoffset,
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC adm data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Read_Reals_at(mbptr.data(), mbcnt, myoffset,
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC adm data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    Kokkos::deep_copy(Kokkos::subview(padm->u_adm, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
    offset_myrank += nout1*nout2*nout3*nadm*sizeof(Real);   // adm u_adm
    myoffset = offset_myrank;
  }

  // read the general-EOS temperature cache, written last (see the size check above),
  // and the derived (p, Gamma_1) cache that follows it
  // THE WARM-START HISTORY comes behind them, and is wanted only by a run that will
  // actually keep one: problem/rt_impl_warm > 0 under the mode-3 column solve.  It is
  // read here but cannot be delivered here -- the Views it belongs in are allocated by
  // the first RT call, long after this -- so it is parked in the staging of
  // utils/two_stream_warm_rst.hpp and the lazy allocation consumes it.  Asking pin
  // rather than two_stream_rt::rt_impl_warm because the problem generator has not run
  // yet and has therefore not set it.
  int warm_want = 0;
  if (pin->DoesParameterExist("problem", "rt_impl_warm") &&
      pin->DoesParameterExist("problem", "rt_implicit_column")) {
    if (pin->GetInteger("problem", "rt_implicit_column") == 3) {
      warm_want = pin->GetInteger("problem", "rt_impl_warm");
    }
  }
  int nwarm_read = (warm_want > 0) ? nwarm_file : 0;
  if (warm_want > 0 && nwarm_file == 0 && global_variable::my_rank == 0) {
    std::cout << "### WARNING: restart file has no mode-3 warm-start history (written "
              << "before it was added, or by a run without one); the first column solve "
              << "cold starts and this restart is not bitwise." << std::endl;
  }
  if (warm_want > 1 && nwarm_file == 1 && global_variable::my_rank == 0) {
    std::cout << "### WARNING: restart file carries only one mode-3 warm-start history "
              << "level and problem/rt_impl_warm = 2 wants two; the first column solve "
              << "does not extrapolate and this restart is not bitwise." << std::endl;
  }
  // the <rad_m1> predictor increment: delivered straight into ipred (allocated by the
  // RadiationM1 constructor) when this run keeps one with the same channel count.  A
  // file without it leaves pred_ok = false: the first step cold starts, as before.
  const bool pred_want = (pradm1 != nullptr) && pradm1->impl_pred;
  const bool pred_read = pred_want && (npred_file > 0) &&
                         (static_cast<int>(pradm1->ipred.extent(1)) == npred_file) &&
                         (static_cast<int>(pradm1->ipred.extent(0)) >= nmb);
  if (pred_want && !pred_read && global_variable::my_rank == 0) {
    std::cout << "### WARNING: restart file has no <rad_m1> implicit_predictor "
              << "increment (written before it was added, or by a run without one); "
              << "the first implicit step cold starts and this restart is not bitwise."
              << std::endl;
  }
  // time_scheme = hesdirk2: the stored slope, when this run wants one of that size.  A
  // file without it leaves t2_ok = false: the first step is backward Euler.
  const bool t2_want = (pradm1 != nullptr) &&
                       (pradm1->time_scheme == radm1::M1_TIME_HESDIRK2);
  const bool t2_read = t2_want && (nt2_file > 0) &&
                       (pradm1->Time2RstNchWant() == nt2_file) &&
                       (static_cast<int>(pradm1->t2k1.extent(0)) >= nmb);
  if (t2_want && !t2_read && global_variable::my_rank == 0) {
    std::cout << "### WARNING: restart file has no <rad_m1> hesdirk2 slope; the first "
              << "step is backward Euler and this restart is not bitwise." << std::endl;
  }
  if (wt_hyd || wt_mhd || nwarm_read > 0 || pred_read || t2_read) {
    const IOWrapperSizeT tail0 = offset_myrank;
    HostArray4D<Real> wtin("rst-wt-in", 1, 1, 1, 1);
    Kokkos::realloc(wtin, nmb, nout3, nout2, nout1);
    // fills wtin with the next per-MeshBlock slab and advances the offsets; the caller
    // copies it where it belongs, which is a 4D array for wtemp and one channel of the
    // 5D wder for the derived cache
    auto read_slab = [&](const char *what) {
      for (int m=0;  m<noutmbs_max; ++m) {
        if (m < noutmbs_min) {
          auto mbptr = Kokkos::subview(wtin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
          int mbcnt = mbptr.size();
          if (resfile.Read_Reals_at_all(mbptr.data(), mbcnt, myoffset,
                                        single_file_per_rank) != mbcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << what << " cache not read correctly from rst file, "
                      << "restart file is broken." << std::endl;
            exit(EXIT_FAILURE);
          }
          myoffset += data_size;
        } else if (m < pm->nmb_thisrank) {
          auto mbptr = Kokkos::subview(wtin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
          int mbcnt = mbptr.size();
          if (resfile.Read_Reals_at(mbptr.data(), mbcnt, myoffset,
                                    single_file_per_rank) != mbcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << what << " cache not read correctly from rst file, "
                      << "restart file is broken." << std::endl;
            exit(EXIT_FAILURE);
          }
          myoffset += data_size;
        }
      }
      offset_myrank += nout1*nout2*nout3*sizeof(Real);
      myoffset = offset_myrank;
    };
    if (wt_hyd) {
      read_slab("hydro");
      Kokkos::deep_copy(Kokkos::subview(phydro->wtemp, std::make_pair(0,nmb), Kokkos::ALL,
                        Kokkos::ALL, Kokkos::ALL), wtin);
    }
    if (wt_mhd) {
      read_slab("mhd");
      Kokkos::deep_copy(Kokkos::subview(pmhd->wtemp, std::make_pair(0,nmb), Kokkos::ALL,
                        Kokkos::ALL, Kokkos::ALL), wtin);
    }
    if (wd_hyd) {
      read_slab("hydro pressure");
      Kokkos::deep_copy(Kokkos::subview(phydro->wder, std::make_pair(0,nmb),
                        static_cast<int>(IDPR), Kokkos::ALL, Kokkos::ALL,
                        Kokkos::ALL), wtin);
      read_slab("hydro Gamma_1");
      Kokkos::deep_copy(Kokkos::subview(phydro->wder, std::make_pair(0,nmb),
                        static_cast<int>(IDG1), Kokkos::ALL, Kokkos::ALL,
                        Kokkos::ALL), wtin);
      // both caches are in place: the first conversion to primitives of this run must
      // not re-derive them, or the restart is not a bitwise continuation.  See
      // Hydro::c2p_freeze_derived.
      phydro->c2p_freeze_derived = true;
    }
    if (wd_mhd) {
      read_slab("mhd pressure");
      Kokkos::deep_copy(Kokkos::subview(pmhd->wder, std::make_pair(0,nmb),
                        static_cast<int>(IDPR), Kokkos::ALL, Kokkos::ALL,
                        Kokkos::ALL), wtin);
      read_slab("mhd Gamma_1");
      Kokkos::deep_copy(Kokkos::subview(pmhd->wder, std::make_pair(0,nmb),
                        static_cast<int>(IDG1), Kokkos::ALL, Kokkos::ALL,
                        Kokkos::ALL), wtin);
      // the MHD counterpart of the line above; see MHD::c2p_freeze_derived
      pmhd->c2p_freeze_derived = true;
    }
    if (nwarm_read > 0) {
      namespace ts = two_stream_rt;
      ts::rt_warm_stage_nmb = nmb;
      ts::rt_warm_stage_n3 = nout3;
      ts::rt_warm_stage_n2 = nout2;
      ts::rt_warm_stage_n1 = nout1;
      ts::rt_warm_stage_nlev = nwarm_read;
      ts::rt_c3_bdt_prev = warm_bdt_prev;
      const std::size_t nw = static_cast<std::size_t>(nmb)*nout3*nout2*nout1;
      read_slab("rt warm start");
      ts::rt_warm_stage.assign(wtin.data(), wtin.data() + nw);
      if (nwarm_read > 1) {
        read_slab("rt warm start 2");
        ts::rt_warm_stage2.assign(wtin.data(), wtin.data() + nw);
      }
    }
    if (pred_read) {
      // step over the warm-start slabs this run did not want
      offset_myrank += (nwarm_file - nwarm_read)*nout1*nout2*nout3*sizeof(Real);
      myoffset = offset_myrank;
      for (int n=0; n<npred_file; ++n) {
        read_slab("rad_m1 predictor");
        Kokkos::deep_copy(Kokkos::subview(pradm1->ipred, std::make_pair(0,nmb), n,
                          Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), wtin);
      }
      pradm1->pred_ok = true;
      pradm1->pred_dt = pred_dt_file;
    }
    if (t2_read) {
      // the slope slabs sit behind every other tail slab of the record
      int nprev = (wt_hyd ? 1 : 0) + (wt_mhd ? 1 : 0) + (wd_hyd ? 2 : 0)
                  + (wd_mhd ? 2 : 0)
                  + nwarm_file + npred_file;
      offset_myrank = tail0 + nprev*nout1*nout2*nout3*sizeof(Real);
      myoffset = offset_myrank;
      for (int n=0; n<nt2_file; ++n) {
        read_slab("rad_m1 hesdirk2 slope");
        pradm1->Time2RstSet(n, wtin, nmb);
      }
      pradm1->t2_ok = true;
      pradm1->pred2_ok = (t2_hv[0] > 0.0);
      pradm1->pred2_dt = t2_hv[0];
      pradm1->t2_dtprev = t2_hv[1];
      pradm1->t2_vprev = (t2_hv[2] > 0.5);
    }
  }

  // call problem generator again to re-initialize data, fn ptrs, as needed
  // second argument true since this IS a restart
  CallProblemGenerator(pin, true);

  // Check that user defined BCs were enrolled if needed
  if (user_bcs) {
    if (user_bcs_func == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "User BCs specified in <mesh> block, but not enrolled "
                << "during restart by SetProblemData()." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  // Check that user defined srcterms were enrolled if needed
  if (user_srcs) {
    if (user_srcs_func == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "User SRCs specified in <problem> block, but not "
                << "enrolled by UserProblem()." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  // Check that user defined history outputs were enrolled if needed
  if (user_hist) {
    if (user_hist_func == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "User history output specified in <problem> block, "
                << "but not enrolled by UserProblem()." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::OutputErrors()
//! \brief Generic function for computing the L1 and L-infty difference between solutions
//! stored in the u0 and u1 registers, and outputting them to an error file.  This is
//! used for linear wave convergence tests, for example.
//! Function requires appropriate solutions already stored in u0 and u1.

void ProblemGenerator::OutputErrors(ParameterInput *pin, Mesh *pm) {
  Real l1_err[16];
  Real linfty_err=0.0;
  int nvars=0,nprev=0;

  // capture class variables for kernel
  auto &indcs = pm->mb_indcs;
  int &nx1 = indcs.nx1;
  int &nx2 = indcs.nx2;
  int &nx3 = indcs.nx3;
  int &is = indcs.is;
  int &js = indcs.js;
  int &ks = indcs.ks;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &size = pmbp->pmb->mb_size;

  // compute errors for Hydro  -----------------------------------------------------------
  if (pmbp->phydro != nullptr) {
    nvars = pmbp->phydro->nhydro;

    auto &is_ideal_ = pmbp->phydro->peos->eos_data.is_ideal;
    auto &u0_ = pmbp->phydro->u0;
    auto &u1_ = pmbp->phydro->u1;

    const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
    const int nkji = nx3*nx2*nx1;
    const int nji  = nx2*nx1;
    array_sum::GlobalSum sum_this_mb;
    Kokkos::parallel_reduce("L1-err",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum, Real &max_err) {
      // compute n,k,j,i indices of thread
      int m = (idx)/nkji;
      int k = (idx - m*nkji)/nji;
      int j = (idx - m*nkji - k*nji)/nx1;
      int i = (idx - m*nkji - k*nji - j*nx1) + is;
      k += ks;
      j += js;

      Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;

      // conserved variables:
      array_sum::GlobalSum evars;
      evars.the_array[IDN] = vol*fabs(u0_(m,IDN,k,j,i) - u1_(m,IDN,k,j,i));
      max_err = fmax(max_err, evars.the_array[IDN]);
      evars.the_array[IM1] = vol*fabs(u0_(m,IM1,k,j,i) - u1_(m,IM1,k,j,i));
      max_err = fmax(max_err, evars.the_array[IM1]);
      evars.the_array[IM2] = vol*fabs(u0_(m,IM2,k,j,i) - u1_(m,IM2,k,j,i));
      max_err = fmax(max_err, evars.the_array[IM2]);
      evars.the_array[IM3] = vol*fabs(u0_(m,IM3,k,j,i) - u1_(m,IM3,k,j,i));
      max_err = fmax(max_err, evars.the_array[IM3]);
      if (is_ideal_) {
        evars.the_array[IEN] = vol*fabs(u0_(m,IEN,k,j,i) - u1_(m,IEN,k,j,i));
        max_err = fmax(max_err, evars.the_array[IEN]);
      }

      // fill rest of the_array with zeros, if narray < NREDUCTION_VARIABLES
      for (int n=nvars; n<NREDUCTION_VARIABLES; ++n) {
        evars.the_array[n] = 0.0;
      }

      // sum into parallel reduce
      mb_sum += evars;
    }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb), Kokkos::Max<Real>(linfty_err));

    // store data into l1_err array
    for (int n=0; n<nvars; ++n) {
      l1_err[n] = sum_this_mb.the_array[n];
    }
    nprev += nvars;
  }

  // compute errors for MHD  -------------------------------------------------------------
  if (pmbp->pmhd != nullptr) {
    nvars = pmbp->pmhd->nmhd + 3;  // include 3-compts of cell-centered B in errors
    auto &is_ideal_ = pmbp->pmhd->peos->eos_data.is_ideal;

    int bindx;
    if (is_ideal_) {
      bindx = 5;
    } else {
      bindx = 4;
    }

    auto &u0_ = pmbp->pmhd->u0;
    auto &u1_ = pmbp->pmhd->u1;
    auto &b0_ = pmbp->pmhd->b0;
    auto &b1_ = pmbp->pmhd->b1;

    const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
    const int nkji = nx3*nx2*nx1;
    const int nji  = nx2*nx1;
    array_sum::GlobalSum sum_this_mb;
    Kokkos::parallel_reduce("L1-err-Sums",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum, Real &max_err) {
      // compute n,k,j,i indices of thread
      int m = (idx)/nkji;
      int k = (idx - m*nkji)/nji;
      int j = (idx - m*nkji - k*nji)/nx1;
      int i = (idx - m*nkji - k*nji - j*nx1) + is;
      k += ks;
      j += js;

      Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;

      // conserved variables:
      array_sum::GlobalSum evars;
      evars.the_array[IDN] = vol*fabs(u0_(m,IDN,k,j,i) - u1_(m,IDN,k,j,i));
      max_err = fmax(max_err, evars.the_array[IDN]);
      evars.the_array[IM1] = vol*fabs(u0_(m,IM1,k,j,i) - u1_(m,IM1,k,j,i));
      max_err = fmax(max_err, evars.the_array[IM1]);
      evars.the_array[IM2] = vol*fabs(u0_(m,IM2,k,j,i) - u1_(m,IM2,k,j,i));
      max_err = fmax(max_err, evars.the_array[IM2]);
      evars.the_array[IM3] = vol*fabs(u0_(m,IM3,k,j,i) - u1_(m,IM3,k,j,i));
      max_err = fmax(max_err, evars.the_array[IM3]);
      if (is_ideal_) {
        evars.the_array[IEN] = vol*fabs(u0_(m,IEN,k,j,i) - u1_(m,IEN,k,j,i));
        max_err = fmax(max_err, evars.the_array[IEN]);
      }

      // cell-centered B
      Real bcc0 = 0.5*(b0_.x1f(m,k,j,i) + b0_.x1f(m,k,j,i+1));
      Real bcc1 = 0.5*(b1_.x1f(m,k,j,i) + b1_.x1f(m,k,j,i+1));
      evars.the_array[bindx] = vol*fabs(bcc0 - bcc1);
      max_err = fmax(max_err, evars.the_array[IEN+1]);

      bcc0 = 0.5*(b0_.x2f(m,k,j,i) + b0_.x2f(m,k,j+1,i));
      bcc1 = 0.5*(b1_.x2f(m,k,j,i) + b1_.x2f(m,k,j+1,i));
      evars.the_array[bindx+1] = vol*fabs(bcc0 - bcc1);
      max_err = fmax(max_err, evars.the_array[IEN+2]);

      bcc0 = 0.5*(b0_.x3f(m,k,j,i) + b0_.x3f(m,k+1,j,i));
      bcc1 = 0.5*(b1_.x3f(m,k,j,i) + b1_.x3f(m,k+1,j,i));
      evars.the_array[bindx+2] = vol*fabs(bcc0 - bcc1);
      max_err = fmax(max_err, evars.the_array[IEN+3]);

      // fill rest of the_array with zeros, if narray < NREDUCTION_VARIABLES
      for (int n=nvars; n<NREDUCTION_VARIABLES; ++n) {
        evars.the_array[n] = 0.0;
      }

      // sum into parallel reduce
      mb_sum += evars;
    }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb), Kokkos::Max<Real>(linfty_err));

    // store data into l1_err array
    for (int n=0; n<nvars; ++n) {
      l1_err[n+nprev] = sum_this_mb.the_array[n];
    }
    nprev += nvars;
  }

#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &l1_err, nprev, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &linfty_err, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif

  // normalize errors by number of cells
  Real vol=  (pmbp->pmesh->mesh_size.x1max - pmbp->pmesh->mesh_size.x1min)
            *(pmbp->pmesh->mesh_size.x2max - pmbp->pmesh->mesh_size.x2min)
            *(pmbp->pmesh->mesh_size.x3max - pmbp->pmesh->mesh_size.x3min);
  for (int i=0; i<nprev; ++i) l1_err[i] = l1_err[i]/vol;
  linfty_err /= vol;

  // compute rms error
  Real rms_err = 0.0;
  for (int i=0; i<nprev; ++i) {
    rms_err += SQR(l1_err[i]);
  }
  rms_err = std::sqrt(rms_err);

  // root process opens output file and writes out errors
  if (global_variable::my_rank == 0) {
    std::string fname;
    fname.assign(pin->GetString("job","basename"));
    fname.append("-errs.dat");
    FILE *pfile;

    // The file exists -- reopen the file in append mode
    if ((pfile = std::fopen(fname.c_str(), "r")) != nullptr) {
      if ((pfile = std::freopen(fname.c_str(), "a", pfile)) == nullptr) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Error output file could not be opened" <<std::endl;
        std::exit(EXIT_FAILURE);
      }

    // The file does not exist -- open the file in write mode and add headers
    } else {
      if ((pfile = std::fopen(fname.c_str(), "w")) == nullptr) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Error output file could not be opened" <<std::endl;
        std::exit(EXIT_FAILURE);
      }
      std::fprintf(pfile, "# Nx1  Nx2  Nx3   Ncycle   RMS-L1       L-infty       ");
      if (pmbp->phydro != nullptr) {
        std::fprintf(pfile,"d_L1          M1_L1         M2_L1         M3_L1         ");
        if (pmbp->phydro->peos->eos_data.is_ideal) {
          std::fprintf(pfile,"E_L1          ");
        }
      }
      if (pmbp->pmhd != nullptr) {
        std::fprintf(pfile,"d_L1          M1_L1         M2_L1         M3_L1         ");
        if (pmbp->pmhd->peos->eos_data.is_ideal) {
          std::fprintf(pfile,"E_L1          ");
        }
        std::fprintf(pfile,"B1_L1         B2_L1         B3_L1");
      }
      std::fprintf(pfile, "\n");
    }

    // write errors
    std::fprintf(pfile, "%04d", pmbp->pmesh->mesh_indcs.nx1);
    std::fprintf(pfile, "  %04d", pmbp->pmesh->mesh_indcs.nx2);
    std::fprintf(pfile, "  %04d", pmbp->pmesh->mesh_indcs.nx3);
    std::fprintf(pfile, "  %05d  %e %e", pmbp->pmesh->ncycle, rms_err, linfty_err);
    for (int i=0; i<nprev; ++i) {
      std::fprintf(pfile, "  %e", l1_err[i]);
    }
    std::fprintf(pfile, "\n");
    std::fclose(pfile);
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::CallProblemGenerator()
//! \brief selects one of the default problem generators compiled automatically with
//! the source code depending on input string in <problem> block ELSE selects a
//! user-defined problem generator function compiled with the code.

void ProblemGenerator::CallProblemGenerator(ParameterInput *pin, bool is_restart) {
#if USER_PROBLEM_ENABLED
  // call user-defined problem generator (if USER_PROBLEM_ENABLED macro defined at build)
  UserProblem(pin, is_restart);
#else
  // else read name of built-in pgen from <problem> block in input file, and call
  std::string pgen_fun_name = pin->GetOrAddString("problem", "pgen_name", "none");

  if (pgen_fun_name.compare("advection") == 0) {
    Advection(pin, is_restart);
  } else if (pgen_fun_name.compare("cpaw") == 0) {
    AlfvenWave(pin, is_restart);
  } else if (pgen_fun_name.compare("gr_bondi") == 0) {
    BondiAccretion(pin, is_restart);
  } else if (pgen_fun_name.compare("cshock") == 0) {
    CShock(pin, is_restart);
  } else if (pgen_fun_name.compare("linear_wave") == 0) {
    LinearWave(pin, is_restart);
  } else if (pgen_fun_name.compare("implode") == 0) {
    LWImplode(pin, is_restart);
  } else if (pgen_fun_name.compare("gr_monopole") == 0) {
    Monopole(pin, is_restart);
  } else if (pgen_fun_name.compare("mri3d") == 0) {
    MRI3d(pin, is_restart);
  } else if (pgen_fun_name.compare("gresho_mhd") == 0) {
    GreshoMHD(pin, is_restart);
  } else if (pgen_fun_name.compare("balsara_vortex") == 0) {
    BalsaraVortex(pin, is_restart);
  } else if (pgen_fun_name.compare("orszag_tang") == 0) {
    OrszagTang(pin, is_restart);
  } else if (pgen_fun_name.compare("rad_linear_wave") == 0) {
    RadiationLinearWave(pin, is_restart);
  } else if (pgen_fun_name.compare("rad_beam") == 0) {
    RadiationBeam(pin, is_restart);
  } else if (pgen_fun_name.compare("rad_m1_beam") == 0) {
    RadiationM1Tests(pin, is_restart);
  } else if (pgen_fun_name.compare("rad_diff2d") == 0) {
    RadDiff2D(pin, is_restart);
  } else if (pgen_fun_name.compare("shock_tube") == 0) {
    ShockTube(pin, is_restart);
  } else if (pgen_fun_name.compare("shwave") == 0) {
    Shwave(pin, is_restart);
  } else if (pgen_fun_name.compare("wb_atm") == 0) {
    WbAtmosphere(pin, is_restart);
  } else if (pgen_fun_name.compare("z4c_boosted_puncture") == 0) {
    Z4cBoostedPuncture(pin, is_restart);
  } else if (pgen_fun_name.compare("z4c_linear_wave") == 0) {
    Z4cLinearWave(pin, is_restart);
  } else if (pgen_fun_name.compare("spherical_collapse") == 0) {
    SphericalCollapse(pin, is_restart);
  } else if (pgen_fun_name.compare("diffusion") == 0) {
    Diffusion(pin, is_restart);
  } else if (pgen_fun_name.compare("gravity") == 0) {
    SelfGravity(pin, is_restart);
  } else if (pgen_fun_name.compare("binary_gravity") == 0) {
    BinaryGravity(pin, is_restart);
  } else if (pgen_fun_name.compare("be_collapse") == 0) {
    BECollapse(pin, is_restart);
  // else, name not set on command line or input file, print warning and quit
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
        << "Problem generator name could not be found in <problem> block in input file"
        << std::endl
        << "and it was not set by -D PROBLEM option on cmake command line during build"
        << std::endl
        << "Rerun cmake with -D PROBLEM=file to specify custom problem generator file"
        << std::endl;;
    std::exit(EXIT_FAILURE);
  }
#endif
  return;
}
