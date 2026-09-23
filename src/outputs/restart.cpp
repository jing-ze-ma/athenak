//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file restart.cpp
//! \brief writes restart files

#include <sys/stat.h>  // mkdir

#include <algorithm>
#include <cstdint>
#include <cstdio>      // fwrite(), fclose(), fopen(), fnprintf(), snprintf()
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility> // make_pair
#include <vector>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "geodesic-grid/geodesic_grid.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "coordinates/adm.hpp"
#include "z4c/compact_object_tracker.hpp"
#include "z4c/z4c.hpp"
#include "radiation/radiation.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"
#include "srcterms/turb_driver.hpp"
#include "pgen/pgen.hpp"
#include "utils/two_stream_warm_rst.hpp"
//#include "outputs.hpp"

//----------------------------------------------------------------------------------------
// constructor: also calls BaseTypeOutput base class constructor

RestartOutput::RestartOutput(ParameterInput *pin, Mesh *pm, OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  // create directories for outputs. Comments in binary.cpp constructor explain why
  mkdir("rst",0775);
  bool single_file_per_rank = op.single_file_per_rank;
  if (single_file_per_rank) {
    char rank_dir[20];
    std::snprintf(rank_dir, sizeof(rank_dir), "rst/rank_%08d/", global_variable::my_rank);
    mkdir(rank_dir, 0775);
  }
}

//----------------------------------------------------------------------------------------
// RestartOutput::LoadOutputData()
// overload of standard load data function specific to restarts.  Loads dependent
// variables, including ghost zones.

void RestartOutput::LoadOutputData(Mesh *pm) {
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
  radiation::Radiation* prad = pm->pmb_pack->prad;
  // <rad_m1>: the four evolved moments (E, F_i), written and read as one more
  // cell-centred block.  Guarded by the pointer everywhere, so a file written by a run
  // WITHOUT <rad_m1> has exactly the layout and the bytes it had before (design sect.6)
  radm1::RadiationM1* pradm1 = pm->pmb_pack->pradm1;
  TurbulenceDriver* pturb=pm->pmb_pack->pturb;
  int nhydro=0, nmhd=0, nrad=0, nm1=0, nforce=3, nadm=0, nz4c=0;
  if (pradm1 != nullptr) {
    nm1 = radm1::M1_NVAR;
  }
  if (phydro != nullptr) {
    nhydro = phydro->nhydro + phydro->nscalars;
  }
  if (pmhd != nullptr) {
    nmhd = pmhd->nmhd + pmhd->nscalars;
  }
  if (pz4c != nullptr) {
    nz4c = pz4c->nz4c;
  } else if (padm != nullptr) {
    nadm = padm->nadm;
  }
  // if the spacetime is evolved, we do not need to checkpoint/recover the ADM variables
  if (prad != nullptr) {
    nrad = prad->prgeo->nangles;
  }

  // Note for restarts, outarrays are dimensioned (m,n,k,j,i)
  if (phydro != nullptr) {
    Kokkos::realloc(outarray_hyd, nmb, nhydro, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_hyd, Kokkos::subview(phydro->u0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (pmhd != nullptr) {
    Kokkos::realloc(outarray_mhd, nmb, nmhd, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_mhd, Kokkos::subview(pmhd->u0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    Kokkos::realloc(outfield.x1f, nmb, nout3, nout2, nout1+1);
    Kokkos::deep_copy(outfield.x1f, Kokkos::subview(pmhd->b0.x1f, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    Kokkos::realloc(outfield.x2f, nmb, nout3, nout2+1, nout1);
    Kokkos::deep_copy(outfield.x2f, Kokkos::subview(pmhd->b0.x2f, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    Kokkos::realloc(outfield.x3f, nmb, nout3+1, nout2, nout1);
    Kokkos::deep_copy(outfield.x3f, Kokkos::subview(pmhd->b0.x3f, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  // the general-EOS temperature cache, see the note on outarray_wth in outputs.hpp
  if (phydro != nullptr && phydro->peos->eos_data.IsGeneral()) {
    Kokkos::realloc(outarray_wth, nmb, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_wth, Kokkos::subview(phydro->wtemp,
                      std::make_pair(0,nmb), Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (pmhd != nullptr && pmhd->peos->eos_data.IsGeneral()) {
    Kokkos::realloc(outarray_wtm, nmb, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_wtm, Kokkos::subview(pmhd->wtemp,
                      std::make_pair(0,nmb), Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  // the general-EOS derived cache p and Gamma_1, see the note on outarray_wdp
  if (phydro != nullptr && phydro->peos->eos_data.IsGeneral()) {
    Kokkos::realloc(outarray_wdp, nmb, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_wdp, Kokkos::subview(phydro->wder,
                      std::make_pair(0,nmb), static_cast<int>(IDPR),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    Kokkos::realloc(outarray_wdg, nmb, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_wdg, Kokkos::subview(phydro->wder,
                      std::make_pair(0,nmb), static_cast<int>(IDG1),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  // and the MHD one, which is the same cache of the same two channels: see the note on
  // outarray_wdpm.  The general-EOS MHD runs (the deep hot Jupiter) are chained restarts
  // and every link boundary was a discontinuity without it.
  if (pmhd != nullptr && pmhd->peos->eos_data.IsGeneral()) {
    Kokkos::realloc(outarray_wdpm, nmb, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_wdpm, Kokkos::subview(pmhd->wder,
                      std::make_pair(0,nmb), static_cast<int>(IDPR),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    Kokkos::realloc(outarray_wdgm, nmb, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_wdgm, Kokkos::subview(pmhd->wder,
                      std::make_pair(0,nmb), static_cast<int>(IDG1),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  // the mode-3 Newton warm-start history, see the note on outarray_wm1
  {
    const int nwm = two_stream_rt::RtWarmLevels();
    if (nwm > 0) {
      Kokkos::realloc(outarray_wm1, nmb, nout3, nout2, nout1);
      Kokkos::deep_copy(outarray_wm1, Kokkos::subview(*two_stream_rt::rt_c3bp_ptr,
                        std::make_pair(0,nmb), Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    }
    if (nwm > 1) {
      Kokkos::realloc(outarray_wm2, nmb, nout3, nout2, nout1);
      Kokkos::deep_copy(outarray_wm2, Kokkos::subview(*two_stream_rt::rt_c3bp2_ptr,
                        std::make_pair(0,nmb), Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    }
  }
  if (prad != nullptr) {
    Kokkos::realloc(outarray_rad, nmb, nrad, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_rad, Kokkos::subview(prad->i0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (pradm1 != nullptr) {
    Kokkos::realloc(outarray_m1, nmb, nm1, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_m1, Kokkos::subview(pradm1->u0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    // milestone 3a: the persistent x1 face fluxes of the implicit scheme
    if (pradm1->transport >= radm1::M1_TRANSPORT_IMPLICIT_X1) {
      Kokkos::realloc(outarray_m1f, nmb, nout3, nout2, nout1+1);
      Kokkos::deep_copy(outarray_m1f, Kokkos::subview(pradm1->f0x1,
                        std::make_pair(0,nmb), Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    }
    // milestone 3b phase B: the transverse face fluxes, in the SAME order the reader
    // (src/pgen/pgen.cpp) takes them back: f0x1, then f0x2, then f0x3.
    if (pradm1->trans_on) {
      Kokkos::realloc(outarray_m1f2, nmb, nout3, nout2+1, nout1);
      Kokkos::deep_copy(outarray_m1f2, Kokkos::subview(pradm1->f0x2,
                        std::make_pair(0,nmb), Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
      if (pradm1->trans_x3) {
        Kokkos::realloc(outarray_m1f3, nmb, nout3+1, nout2, nout1);
        Kokkos::deep_copy(outarray_m1f3, Kokkos::subview(pradm1->f0x3,
                          std::make_pair(0,nmb), Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
      }
    }
    // implicit_predictor = step: the stored previous-step increment (radm1::ipred)
    if (pradm1->impl_pred && pradm1->pred_ok) {
      const int npc = static_cast<int>(pradm1->ipred.extent(1));
      Kokkos::realloc(outarray_m1p, nmb, npc, nout3, nout2, nout1);
      Kokkos::deep_copy(outarray_m1p, Kokkos::subview(pradm1->ipred,
                        std::make_pair(0,nmb), Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                        Kokkos::ALL));
    }
    // time_scheme = hesdirk2: the FSAL slope and what goes with it
    if (pradm1->Time2RstNch() > 0) {
      const int nct = pradm1->Time2RstNch();
      DvceArray5D<Real> tmp("rst-m1t", nmb, nct, nout3, nout2, nout1);
      pradm1->Time2RstPack(tmp, nmb);
      Kokkos::realloc(outarray_m1t, nmb, nct, nout3, nout2, nout1);
      Kokkos::deep_copy(outarray_m1t, tmp);
    }
  }
  if (pturb != nullptr) {
    Kokkos::realloc(outarray_force, nmb, nforce, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_force, Kokkos::subview(pturb->force, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (pz4c != nullptr) {
    Kokkos::realloc(outarray_z4c, nmb, nz4c, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_z4c, Kokkos::subview(pz4c->u0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  } else if (padm != nullptr) {
    Kokkos::realloc(outarray_adm, nmb, nadm, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_adm, Kokkos::subview(padm->u_adm, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }

  // calculate max/min number of MeshBlocks across all ranks
  noutmbs_max = pm->nmb_eachrank[0];
  noutmbs_min = pm->nmb_eachrank[0];
  for (int i=0; i<(global_variable::nranks); ++i) {
    noutmbs_max = std::max(noutmbs_max,pm->nmb_eachrank[i]);
    noutmbs_min = std::min(noutmbs_min,pm->nmb_eachrank[i]);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RestartOutput:::WriteOutputFile(Mesh *pm)
//  \brief Cycles over all MeshBlocks and writes everything to a single restart file

void RestartOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  // get spatial dimensions of arrays, including ghost zones
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int nout1 = indcs.nx1 + 2*(indcs.ng);
  int nout2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int nout3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  hydro::Hydro* phydro = pm->pmb_pack->phydro;
  mhd::MHD* pmhd = pm->pmb_pack->pmhd;
  radiation::Radiation* prad = pm->pmb_pack->prad;
  radm1::RadiationM1* pradm1 = pm->pmb_pack->pradm1;
  TurbulenceDriver* pturb=pm->pmb_pack->pturb;
  z4c::Z4c* pz4c = pm->pmb_pack->pz4c;
  adm::ADM* padm = pm->pmb_pack->padm;
  // the OPTIONAL, MARKED pgen state block (see ProblemGenerator::pgen_rst_write_func).
  // Empty -- which is every problem generator but red_giant with problem/mlt_alpha > 0 --
  // and not one byte is written, so the file layout is exactly what it always was.
  std::vector<char> pgen_state;
  if (pm->pgen != nullptr && pm->pgen->pgen_rst_write_func != nullptr) {
    pgen_state = (pm->pgen->pgen_rst_write_func)();
  }
  // the mode-3 Newton warm-start history: how many levels this run keeps (0 when the
  // warm start is off, or on but not yet allocated), and the little marked header that
  // says so -- int32 nlev, int32 pad, then the previous call's bdt, which is the
  // denominator of the rt_impl_warm = 2 extrapolation ratio and is state as much as the
  // levels are.  See utils/two_stream_warm_rst.hpp.
  const int nwarm = two_stream_rt::RtWarmLevels();
  char warm_hdr[2*sizeof(std::int32_t) + sizeof(Real)];
  {
    const std::int32_t hdr[2] = {static_cast<std::int32_t>(nwarm), 0};
    const Real bdtp = two_stream_rt::rt_c3_bdt_prev;
    std::memcpy(&(warm_hdr[0]), &(hdr[0]), sizeof(hdr));
    std::memcpy(&(warm_hdr[0]) + sizeof(hdr), &bdtp, sizeof(bdtp));
  }
  // <rad_m1> implicit_predictor = step: the number of ipred channels written (0 when no
  // increment is stored, and then not one byte), and its marked header -- int32 have,
  // int32 ncomp, Real pred_dt.  See radm1::kM1PredRstMagic.
  const int npred = (pradm1 != nullptr && pradm1->impl_pred && pradm1->pred_ok) ?
                    static_cast<int>(pradm1->ipred.extent(1)) : 0;
  char pred_hdr[2*sizeof(std::int32_t) + sizeof(Real)];
  {
    const std::int32_t hdr[2] = {(npred > 0) ? 1 : 0, static_cast<std::int32_t>(npred)};
    const Real pdt = (pradm1 != nullptr) ? pradm1->pred_dt : 0.0;
    std::memcpy(&(pred_hdr[0]), &(hdr[0]), sizeof(hdr));
    std::memcpy(&(pred_hdr[0]) + sizeof(hdr), &pdt, sizeof(pdt));
  }
  // <rad_m1> time_scheme = hesdirk2: the slope block, behind the predictor one --
  // int32 have, int32 nch, Real pred2_dt, Real dt_prev, Real vprev
  const int nt2 = (pradm1 != nullptr) ? pradm1->Time2RstNch() : 0;
  // <rad_m1> implicit_one_pass: its state, behind the hesdirk2 header
  const bool wonep = (pradm1 != nullptr) && (pradm1->impl_onep > 0);
  Real onep_hdr[9];
  if (wonep) {
    for (int t = 0; t < 3; ++t) {
      onep_hdr[t] = pradm1->onep_qa[t];
      onep_hdr[3+t] = pradm1->onep_qb[t];
      onep_hdr[6+t] = pradm1->onep_cnt[t];
    }
  }
  char t2_hdr[2*sizeof(std::int32_t) + 3*sizeof(Real)];
  {
    const std::int32_t hdr[2] = {(nt2 > 0) ? 1 : 0, static_cast<std::int32_t>(nt2)};
    Real hv[3] = {0.0, 0.0, 0.0};
    if (pradm1 != nullptr) {
      hv[0] = pradm1->pred2_ok ? pradm1->pred2_dt : 0.0;
      hv[1] = pradm1->t2_dtprev;
      hv[2] = pradm1->t2_vprev ? 1.0 : 0.0;
    }
    std::memcpy(&(t2_hdr[0]), &(hdr[0]), sizeof(hdr));
    std::memcpy(&(t2_hdr[0]) + sizeof(hdr), &(hv[0]), sizeof(hv));
  }
  int nhydro=0, nmhd=0, nrad=0, nm1=0, nforce=3, nz4c=0, nadm=0, nco=0;
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
    nco = pz4c->ptracker.size();
  } else if (padm != nullptr) {
    nadm = padm->nadm;
  }
  bool single_file_per_rank = out_params.single_file_per_rank;
  std::string fname;
  if (single_file_per_rank) {
    // Generate a directory and filename for each rank
    // create filename: "rst/rank_YYYYYYY/file_basename" + "." + XXXXX + ".rst"
    // where YYYYYYY = 8-digit rank number
    // where XXXXX = 5-digit file_number
    char rank_dir[20];
    char number[7];
    std::snprintf(number, sizeof(number), ".%05d", out_params.file_number);
    std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d/", global_variable::my_rank);
    fname = std::string("rst/") + std::string(rank_dir) + out_params.file_basename
      + number + ".rst";

    // Debugging output to check directory and filename
    // std::cout << "Rank " << global_variable::my_rank << " generated filename: "
    //           << fname << std::endl;
  } else {
    // Existing behavior: single restart file
    // create filename: "rst/file_basename" + "." + XXXXX + ".rst"
    // where XXXXX = 5-digit file_number
    char number[7];
    std::snprintf(number, sizeof(number), ".%05d", out_params.file_number);
    fname = std::string("rst/") + out_params.file_basename + number + ".rst";
  }
  // increment counters now so values for *next* dump are stored in restart file
  out_params.file_number++;
  if (out_params.last_time < 0.0) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  pin->SetInteger(out_params.block_name, "file_number", out_params.file_number);
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);

  // create string holding input parameters (copy of input file)
  std::stringstream ost;
  pin->ParameterDump(ost);
  std::string sbuf = ost.str();

  //--- STEP 1.  Root process writes header data (input file, critical variables)
  // Input file data is read by ParameterInput on restart, and the remaining header
  // variables are read in Mesh::BuildTreeFromRestart()

  // open file and  write the header; this part is serial
  IOWrapper resfile;
  resfile.Open(fname.c_str(), IOWrapper::FileMode::write, single_file_per_rank);
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    // output the input parameters (input file)
    resfile.Write_any_type(sbuf.c_str(), sbuf.size(), "byte", single_file_per_rank);

    // output Mesh information
    resfile.Write_any_type(&(pm->nmb_total), (sizeof(int)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->root_level), (sizeof(int)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->mesh_size), (sizeof(RegionSize)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->mesh_indcs), (sizeof(RegionIndcs)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->mb_indcs), (sizeof(RegionIndcs)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->time), (sizeof(Real)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->dt), (sizeof(Real)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->ncycle), (sizeof(int)), "byte",
                            single_file_per_rank);
  }
  //--- STEP 2.  Root process writes list of logical locations and cost of MeshBlocks
  // This data read in Mesh::BuildTreeFromRestart()

  if (global_variable::my_rank == 0 || single_file_per_rank) {
    resfile.Write_any_type(&(pm->lloc_eachmb[0]),(pm->nmb_total)*sizeof(LogicalLocation),
                           "byte", single_file_per_rank);
    resfile.Write_any_type(&(pm->cost_eachmb[0]), (pm->nmb_total)*sizeof(float),
                           "byte", single_file_per_rank);
  }

  //--- STEP 3.  Root process writes internal state of objects that require it
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    // store z4c information
    if (pz4c != nullptr) {
      resfile.Write_any_type(&(pz4c->last_output_time), sizeof(Real), "byte",
                             single_file_per_rank);
    }
    // output puncture tracker data
    if (nco > 0) {
      for (auto & pt : pz4c->ptracker) {
        resfile.Write_any_type(pt->GetPos(), 3*sizeof(Real), "byte",
                               single_file_per_rank);
      }
    }
    // turbulence driver internal RNG
    if (pturb != nullptr) {
      resfile.Write_any_type(&(pturb->rstate), sizeof(RNG_State), "byte",
                             single_file_per_rank);
    }
    // the pgen's own relaxed state, marker then length then payload
    if (!pgen_state.empty()) {
      IOWrapperSizeT nb = pgen_state.size();
      resfile.Write_any_type(&(kPgenRstMagic[0]), sizeof(kPgenRstMagic), "byte",
                             single_file_per_rank);
      resfile.Write_any_type(&nb, sizeof(IOWrapperSizeT), "byte", single_file_per_rank);
      resfile.Write_any_type(pgen_state.data(), nb, "byte", single_file_per_rank);
    }
    // the mode-3 warm-start header, marker then length then payload, read by the same
    // eight-byte peek the pgen state block uses (utils/two_stream_warm_rst.hpp)
    if (nwarm > 0) {
      IOWrapperSizeT nb = sizeof(warm_hdr);
      resfile.Write_any_type(&(two_stream_rt::kRtWarmRstMagic[0]),
                             sizeof(two_stream_rt::kRtWarmRstMagic), "byte",
                             single_file_per_rank);
      resfile.Write_any_type(&nb, sizeof(IOWrapperSizeT), "byte", single_file_per_rank);
      resfile.Write_any_type(&(warm_hdr[0]), nb, "byte", single_file_per_rank);
    }
    // the implicit_predictor header, same marked form, behind the warm-start one
    if (npred > 0) {
      IOWrapperSizeT nb = sizeof(pred_hdr);
      resfile.Write_any_type(&(radm1::kM1PredRstMagic[0]),
                             sizeof(radm1::kM1PredRstMagic), "byte",
                             single_file_per_rank);
      resfile.Write_any_type(&nb, sizeof(IOWrapperSizeT), "byte", single_file_per_rank);
      resfile.Write_any_type(&(pred_hdr[0]), nb, "byte", single_file_per_rank);
    }
    // the hesdirk2 slope header, same marked form, behind the predictor one
    if (nt2 > 0) {
      IOWrapperSizeT nb = sizeof(t2_hdr);
      resfile.Write_any_type(&(radm1::kM1Time2RstMagic[0]),
                             sizeof(radm1::kM1Time2RstMagic), "byte",
                             single_file_per_rank);
      resfile.Write_any_type(&nb, sizeof(IOWrapperSizeT), "byte", single_file_per_rank);
      resfile.Write_any_type(&(t2_hdr[0]), nb, "byte", single_file_per_rank);
    }
    if (wonep) {
      IOWrapperSizeT nb = sizeof(onep_hdr);
      resfile.Write_any_type(&(radm1::kM1OnePassRstMagic[0]),
                             sizeof(radm1::kM1OnePassRstMagic), "byte",
                             single_file_per_rank);
      resfile.Write_any_type(&nb, sizeof(IOWrapperSizeT), "byte", single_file_per_rank);
      resfile.Write_any_type(&(onep_hdr[0]), nb, "byte", single_file_per_rank);
    }
  }

  //--- STEP 4.  All ranks write data over all MeshBlocks (5D arrays) in parallel
  // This data read in ProblemGenerator constructor for restarts

  // total size of all cell-centered variables and face-centered fields to be written by
  // this rank
  // the general-EOS temperature cache is appended at the END of each MeshBlock record,
  // so that a file written without it (any run with an ideal gas, and every file written
  // before this was added) differs only in its tail and is still readable -- see the
  // size check in pgen.cpp's restart constructor.
  bool wt_hyd = (phydro != nullptr) && phydro->peos->eos_data.IsGeneral();
  bool wt_mhd = (pmhd != nullptr) && pmhd->peos->eos_data.IsGeneral();
  IOWrapperSizeT data_size = 0;
  if (phydro != nullptr) {
    data_size += nout1*nout2*nout3*nhydro*sizeof(Real); // hydro u0
  }
  if (pmhd != nullptr) {
    data_size += nout1*nout2*nout3*nmhd*sizeof(Real);   // mhd u0
    data_size += (nout1+1)*nout2*nout3*sizeof(Real);    // mhd b0.x1f
    data_size += nout1*(nout2+1)*nout3*sizeof(Real);    // mhd b0.x2f
    data_size += nout1*nout2*(nout3+1)*sizeof(Real);    // mhd b0.x3f
  }
  if (prad != nullptr) {
    data_size += nout1*nout2*nout3*nrad*sizeof(Real);   // radiation i0
  }
  if (pradm1 != nullptr) {
    data_size += nout1*nout2*nout3*nm1*sizeof(Real);    // rad_m1 u0
    if (pradm1->transport >= radm1::M1_TRANSPORT_IMPLICIT_X1) {
      data_size += (nout1+1)*nout2*nout3*sizeof(Real);  // rad_m1 f0x1 (milestone 3a)
    }
    if (pradm1->trans_on) {
      data_size += nout1*(nout2+1)*nout3*sizeof(Real);  // rad_m1 f0x2 (3b phase B)
      if (pradm1->trans_x3) {
        data_size += nout1*nout2*(nout3+1)*sizeof(Real);  // rad_m1 f0x3
      }
    }
  }
  if (pturb != nullptr) {
    data_size += nout1*nout2*nout3*nforce*sizeof(Real); // forcing
  }
  if (pz4c != nullptr) {
    data_size += nout1*nout2*nout3*nz4c*sizeof(Real);   // z4c u0
  } else if (padm != nullptr) {
    data_size += nout1*nout2*nout3*nadm*sizeof(Real);   // adm u_adm
  }
  if (wt_hyd) {
    data_size += nout1*nout2*nout3*sizeof(Real);        // hydro wtemp
  }
  if (wt_mhd) {
    data_size += nout1*nout2*nout3*sizeof(Real);        // mhd wtemp
  }
  if (wt_hyd) {
    data_size += 2*nout1*nout2*nout3*sizeof(Real);      // hydro wder (IDPR, IDG1)
  }
  if (wt_mhd) {
    data_size += 2*nout1*nout2*nout3*sizeof(Real);      // mhd wder (IDPR, IDG1)
  }
  if (nwarm > 0) {
    data_size += nwarm*nout1*nout2*nout3*sizeof(Real);  // rt_c3bp (+ rt_c3bp2)
  }
  if (npred > 0) {
    data_size += npred*nout1*nout2*nout3*sizeof(Real);  // rad_m1 ipred
  }
  if (nt2 > 0) {
    data_size += nt2*nout1*nout2*nout3*sizeof(Real);    // rad_m1 hesdirk2 slope
  }
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    resfile.Write_any_type(&(data_size), sizeof(IOWrapperSizeT), "byte",
                            single_file_per_rank);
  }

  // calculate size of data written in Steps 1-2 above
  IOWrapperSizeT step1size = sbuf.size()*sizeof(char) + 3*sizeof(int) + 2*sizeof(Real) +
                             sizeof(RegionSize) + 2*sizeof(RegionIndcs);
  IOWrapperSizeT step2size = (pm->nmb_total)*(sizeof(LogicalLocation) + sizeof(float));

  IOWrapperSizeT step3size = 3*nco*sizeof(Real);
  if (pz4c != nullptr) step3size += sizeof(Real);
  if (pturb != nullptr) step3size += sizeof(RNG_State);
  if (!pgen_state.empty()) {
    step3size += sizeof(kPgenRstMagic) + sizeof(IOWrapperSizeT) + pgen_state.size();
  }
  if (nwarm > 0) {
    step3size += sizeof(two_stream_rt::kRtWarmRstMagic) + sizeof(IOWrapperSizeT)
                 + sizeof(warm_hdr);
  }
  if (npred > 0) {
    step3size += sizeof(radm1::kM1PredRstMagic) + sizeof(IOWrapperSizeT)
                 + sizeof(pred_hdr);
  }
  if (nt2 > 0) {
    step3size += sizeof(radm1::kM1Time2RstMagic) + sizeof(IOWrapperSizeT)
                 + sizeof(t2_hdr);
  }
  if (wonep) {
    step3size += sizeof(radm1::kM1OnePassRstMagic) + sizeof(IOWrapperSizeT)
                 + sizeof(onep_hdr);
  }

  // write cell-centered variables in parallel
  IOWrapperSizeT offset_myrank = (step1size + step2size + step3size
                                  + sizeof(IOWrapperSizeT));

  if (!single_file_per_rank) {
    offset_myrank += data_size*(pm->gids_eachrank[global_variable::my_rank]);
  }

  IOWrapperSizeT myoffset = offset_myrank;

  // write cell-centered variables, one MeshBlock at a time (but parallelized over all
  // ranks). MeshBlocks are written seperately to reduce number of data elements per write
  // call, to avoid exceeding 2^31 limit for very large grids per MPI rank.
  if (phydro != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_hyd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered hydro data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_hyd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered hydro data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nhydro*sizeof(Real); // hydro u0
    myoffset = offset_myrank;
  }
  if (pmhd != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_mhd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered mhd data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_mhd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered mhd data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nmhd*sizeof(Real);   // mhd u0
    myoffset = offset_myrank;

    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to x1-face field
        auto x1fptr = Kokkos::subview(outfield.x1f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        int fldcnt = x1fptr.size();
        if (resfile.Write_any_type_at_all(x1fptr.data(),fldcnt,myoffset,"Real",
                                          single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x1f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        // get ptr to x2-face field
        auto x2fptr = Kokkos::subview(outfield.x2f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = x2fptr.size();
        if (resfile.Write_any_type_at_all(x2fptr.data(),fldcnt,myoffset,"Real",
                                          single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x2f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        // get ptr to x3-face field
        auto x3fptr = Kokkos::subview(outfield.x3f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = x3fptr.size();
        if (resfile.Write_any_type_at_all(x3fptr.data(),fldcnt,myoffset,"Real",
                                          single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x3f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        myoffset += data_size-(x1fptr.size()+x2fptr.size()+x3fptr.size())*sizeof(Real);

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to x1-face field
        auto x1fptr = Kokkos::subview(outfield.x1f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        int fldcnt = x1fptr.size();
        if (resfile.Write_any_type_at(x1fptr.data(),fldcnt,myoffset,"Real",
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x1f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        // get ptr to x2-face field
        auto x2fptr = Kokkos::subview(outfield.x2f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = x2fptr.size();
        if (resfile.Write_any_type_at(x2fptr.data(),fldcnt,myoffset,"Real",
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x2f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        // get ptr to x3-face field
        auto x3fptr = Kokkos::subview(outfield.x3f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = x3fptr.size();
        if (resfile.Write_any_type_at(x3fptr.data(),fldcnt,myoffset,"Real",
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x3f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        myoffset += data_size-(x1fptr.size()+x2fptr.size()+x3fptr.size())*sizeof(Real);
      }
    }
    offset_myrank += (nout1+1)*nout2*nout3*sizeof(Real);    // mhd b0.x1f
    offset_myrank += nout1*(nout2+1)*nout3*sizeof(Real);    // mhd b0.x2f
    offset_myrank += nout1*nout2*(nout3+1)*sizeof(Real);    // mhd b0.x3f
    myoffset = offset_myrank;
  }

  if (prad != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_rad, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered rad data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_rad, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(),mbcnt,myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered rad data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nrad*sizeof(Real);   // radiation i0
    myoffset = offset_myrank;
  }

  if (pradm1 != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        auto mbptr = Kokkos::subview(outarray_m1, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered rad_m1 data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        auto mbptr = Kokkos::subview(outarray_m1, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(),mbcnt,myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered rad_m1 data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nm1*sizeof(Real);    // rad_m1 u0
    myoffset = offset_myrank;

    // milestone 3a / 3b phase B: the persistent face fluxes, same loop, same order.
    // ONE writer for the three face arrays, so that the order here cannot drift from the
    // order of the reader in src/pgen/pgen.cpp.
    auto wrface = [&](HostArray4D<Real> &arr) {
      for (int m=0;  m<noutmbs_max; ++m) {
        if (m < noutmbs_min) {
          auto mbptr = Kokkos::subview(arr, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
          int mbcnt = mbptr.size();
          if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                            single_file_per_rank) != mbcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "rad_m1 face data not written correctly to rst file, "
            << "restart file is broken." << std::endl;
            exit(EXIT_FAILURE);
          }
          myoffset += data_size;
        } else if (m < pm->nmb_thisrank) {
          auto mbptr = Kokkos::subview(arr, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
          int mbcnt = mbptr.size();
          if (resfile.Write_any_type_at(mbptr.data(),mbcnt,myoffset,"Real",
                                        single_file_per_rank) != mbcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "rad_m1 face data not written correctly"
                      << " to rst file, restart file is broken." << std::endl;
            exit(EXIT_FAILURE);
          }
          myoffset += data_size;
        }
      }
    };
    if (pradm1->transport >= radm1::M1_TRANSPORT_IMPLICIT_X1) {
      for (int m=0;  m<noutmbs_max; ++m) {
        if (m < noutmbs_min) {
          auto mbptr = Kokkos::subview(outarray_m1f, m, Kokkos::ALL, Kokkos::ALL,
                                       Kokkos::ALL);
          int mbcnt = mbptr.size();
          if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                            single_file_per_rank) != mbcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "rad_m1 face data not written correctly to rst file, "
            << "restart file is broken." << std::endl;
            exit(EXIT_FAILURE);
          }
          myoffset += data_size;
        } else if (m < pm->nmb_thisrank) {
          auto mbptr = Kokkos::subview(outarray_m1f, m, Kokkos::ALL, Kokkos::ALL,
                                       Kokkos::ALL);
          int mbcnt = mbptr.size();
          if (resfile.Write_any_type_at(mbptr.data(),mbcnt,myoffset,"Real",
                                        single_file_per_rank) != mbcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "rad_m1 face data not written correctly"
                      << " to rst file, restart file is broken." << std::endl;
            exit(EXIT_FAILURE);
          }
          myoffset += data_size;
        }
      }
      offset_myrank += (nout1+1)*nout2*nout3*sizeof(Real);   // rad_m1 f0x1
      myoffset = offset_myrank;
    }
    if (pradm1->trans_on) {
      wrface(outarray_m1f2);
      offset_myrank += nout1*(nout2+1)*nout3*sizeof(Real);   // rad_m1 f0x2
      myoffset = offset_myrank;
      if (pradm1->trans_x3) {
        wrface(outarray_m1f3);
        offset_myrank += nout1*nout2*(nout3+1)*sizeof(Real);   // rad_m1 f0x3
        myoffset = offset_myrank;
      }
    }
  }

  if (pturb != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_force, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered turb data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_force, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered turb data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nforce*sizeof(Real); // forcing
    myoffset = offset_myrank;
  }

  if (pz4c != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_z4c, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered z4c data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_z4c, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered z4c data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nz4c*sizeof(Real); // z4c u0
    myoffset = offset_myrank;
  } else if (padm != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_adm, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered adm data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_adm, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered adm data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nadm*sizeof(Real); // adm u_adm
    myoffset = offset_myrank;
  }

  // write the general-EOS temperature cache last, one MeshBlock at a time, exactly as
  // the arrays above.  Same loop for hydro and MHD, so it is written once here.
  auto write_wtemp = [&](const auto &a, const char *what) {
    for (int m=0;  m<noutmbs_max; ++m) {
      if (m < noutmbs_min) {
        auto mbptr = Kokkos::subview(a, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << what << " cache not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      } else if (m < pm->nmb_thisrank) {
        auto mbptr = Kokkos::subview(a, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << what << " cache not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*sizeof(Real);
    myoffset = offset_myrank;
  };
  if (wt_hyd) { write_wtemp(outarray_wth, "hydro"); }
  if (wt_mhd) { write_wtemp(outarray_wtm, "mhd"); }
  // the derived cache follows it, same layout, same loop (see outarray_wdp)
  if (wt_hyd) {
    write_wtemp(outarray_wdp, "hydro pressure");
    write_wtemp(outarray_wdg, "hydro Gamma_1");
  }
  if (wt_mhd) {
    write_wtemp(outarray_wdpm, "mhd pressure");
    write_wtemp(outarray_wdgm, "mhd Gamma_1");
  }
  // and the warm-start history behind that, same layout, same loop
  if (nwarm > 0) { write_wtemp(outarray_wm1, "rt warm start"); }
  if (nwarm > 1) { write_wtemp(outarray_wm2, "rt warm start 2"); }
  // and the implicit_predictor increment behind that, one slab per channel
  for (int n=0; n<npred; ++n) {
    write_wtemp(Kokkos::subview(outarray_m1p, Kokkos::ALL, n, Kokkos::ALL, Kokkos::ALL,
                                Kokkos::ALL), "rad_m1 predictor");
  }
  // and the hesdirk2 slope behind that
  for (int n=0; n<nt2; ++n) {
    write_wtemp(Kokkos::subview(outarray_m1t, Kokkos::ALL, n, Kokkos::ALL, Kokkos::ALL,
                                Kokkos::ALL), "rad_m1 hesdirk2 slope");
  }

  // close file, clean up
  resfile.Close(single_file_per_rank);

  // BITWISE RESTARTS: the dynamic well-balanced background (hydro/mhd wbq0) is a CACHE
  // rebuilt only every wb_cache_every cycles, and it is NOT restart state.  A restarted
  // run rebuilds it from the state in this file, while the straight run keeps one built
  // up to wb_cache_every-1 cycles earlier.  Invalidate it here so both runs rebuild it
  // from exactly the state this file holds, on the very next stage 1.
  if (pm->pmb_pack->phydro != nullptr) { pm->pmb_pack->phydro->wb_cache_built = false; }
  if (pm->pmb_pack->pmhd != nullptr) { pm->pmb_pack->pmhd->wb_cache_built = false; }

  return;
}
