//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file build_tree.cpp
//! \brief Functions to build MeshBlockTreee, both for new runs and restarts

#include <iostream>
#include <cinttypes>
#include <limits> // numeric_limits<>
#include <memory> // make_unique<>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh.hpp"
#include "coordinates/cell_locations.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn CheckCubedSphereRefinement
//! \brief Refuse a refined cubed-sphere mesh whose level boundaries touch a panel seam.
//!
//! Refinement is supported on the cubed sphere only where the coarse/fine interfaces are
//! RADIAL. A level boundary that lies on a seam would have to prolongate and restrict
//! ACROSS the seam, and the halo there is not a copy: it applies a tangent-basis
//! transform and an along-seam resample (see bvals_cc.cpp). Neither operator does, so
//! the ghost data would be wrong by O(1) in the velocity and the run would look
//! plausible while being inconsistent. Checked here rather than from the input file
//! because it is a property of the built tree, and refuses rather than produces a
//! quietly wrong grid.

void CheckCubedSphereRefinement(MeshBlockPack *pmbp) {
  auto &nghbr = pmbp->pmb->nghbr;
  auto &mblev = pmbp->pmb->mb_lev;
  auto &mbpanel = pmbp->pmb->mb_panel;
  nghbr.sync_host();
  mblev.sync_host();
  mbpanel.sync_host();
  for (int m=0; m<pmbp->nmb_thispack; ++m) {
    for (int n=0; n<pmbp->pmb->nnghbr; ++n) {
      if (nghbr.h_view(m,n).gid < 0) continue;
      if (nghbr.h_view(m,n).panel == mbpanel.h_view(m)) continue;
      // FACE buffers only. x1 faces (n < 8) are radial and never cross a seam; the x2
      // and x3 face buffers are 8-15 and 24-31. An EDGE or CORNER neighbour can be both
      // across a seam and across a level boundary even when every level boundary is
      // radial -- the block at the outer edge of a refined shell meets, diagonally, the
      // coarse block of the next panel -- and refusing those would rule out radial
      // refinement altogether. Those ng x ng buffers are plain copies that a
      // dimensionally split PLM + HLLC sweep never reconstructs through, which is the
      // same ground on which the halo already leaves them untransformed
      // (see the note in bvals_cc.cpp). A higher-order reconstruction would need them.
      const bool is_face = ((n >= 8 && n < 16) || (n >= 24 && n < 32));
      if (!is_face) continue;
      if (nghbr.h_view(m,n).lev != mblev.h_view(m)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "mesh/use_cubed_sphere supports refinement only where the "
          << "coarse/fine interfaces are RADIAL, but MeshBlock gid="
          << pmbp->pmb->mb_gid.h_view(m) << " on panel " << mbpanel.h_view(m)
          << " (level " << mblev.h_view(m) << ") meets a block on panel "
          << nghbr.h_view(m,n).panel << " at level " << nghbr.h_view(m,n).lev
          << " across a SEAM. Prolongation and restriction have no seam path."
          << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::BuildTreeFromScratch():
//! Constructs MeshBlockTree, creates MeshBlockPack (containing the physics modules), and
//! divides grid into MeshBlock(s) for new runs (starting from scratch), using parameters
//! read from input file.  Also does initial load balance based on simple cost estimate.

void Mesh::BuildTreeFromScratch(ParameterInput *pin) {
  // calculate the number of MeshBlocks at root level in each dir
  nmb_rootx1 = mesh_indcs.nx1/mb_indcs.nx1;
  nmb_rootx2 = mesh_indcs.nx2/mb_indcs.nx2;
  nmb_rootx3 = mesh_indcs.nx3/mb_indcs.nx3;

  if (use_cubed_sphere) {
    // The two PANEL-TANGENTIAL axes are x2 (xi) and x3 (eta) -- x1 is radial -- and a
    // seam glues an x2 face of one panel to an x3 face of another, so the root grid must
    // carry the same number of MeshBlocks along each. The radial direction is free.
    if (nmb_rootx2 != nmb_rootx3) {
      std::cout << "### FATAL ERROR: Cubed sphere requires the same number of MeshBlocks "
                << "along x2 and x3 (the two panel-tangential axes)\n";
      std::exit(EXIT_FAILURE);
    }
  }

  // find maximum number of MeshBlocks at root level in any dir
  int nmbmax = (nmb_rootx1 > nmb_rootx2) ? nmb_rootx1 : nmb_rootx2;
  nmbmax = (nmbmax > nmb_rootx3) ? nmbmax : nmb_rootx3;

  // Find smallest N such that 2^N > max number of MeshBlocks in any dimension (nmbmax)
  // Then N is logical level of root grid.  2^N implemented as left-shift (1<<root_level)
  for (root_level=0; ((1<<root_level) < nmbmax); root_level++) {}
  int current_level = root_level;

  // Construct tree and create root grid
  panel_trees.resize(npanels);
  for (int p = 0; p < npanels; ++p) {
    panel_trees[p] = std::make_unique<MeshBlockTree>(this);
    // store panel id somewhere
    panel_trees[p]->SetPanelID(p);
    panel_trees[p]->CreateRootGrid();
  }
  ptree = panel_trees[0].get();

  // Error check properties of input paraemters for SMR/AMR meshes.
  if (adaptive) {
    // Refinement walks a single tree through ptree, which is panel 0. With six panels
    // that would silently refine one sixth of the mesh and leave the rest alone, so
    // refuse rather than produce a quietly wrong grid.
    if (npanels > 1) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Adaptive mesh refinement is not implemented for the "
                << "cubed sphere: refinement operates on a single MeshBlockTree, but a "
                << "cubed-sphere mesh has one tree per panel." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    max_level = pin->GetOrAddInteger("mesh_refinement", "num_levels", 1) + root_level - 1;
    if (max_level > 31) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Number of refinement levels must be smaller than "
                << 31 - root_level + 1 << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {
    max_level = 31;
  }

  // Read <refined_region> blocks and construct tree accordingly
  // These regions can be used with both SMR (in which case they will remain fixed) and
  // AMR (in which case they may be defined, unless the location refinement criteria used)
  if (multilevel) {
    // error check that number of cells in MeshBlock divisible by two
    if (mb_indcs.nx1 % 2 != 0 ||
       (mb_indcs.nx2 % 2 != 0 && multi_d) ||
       (mb_indcs.nx3 % 2 != 0 && three_d)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Number of cells in MeshBlock must be divisible by 2 "
                << "with SMR or AMR." << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // cycle through ParameterInput list and find "refined_region" blocks, extract data
    // and expand MeshBlockTree
    for (auto it = pin->block.begin(); it != pin->block.end(); ++it) {
      if (it->block_name.compare(0, 14, "refined_region") == 0) {
        RegionSize ref_size;
        ref_size.x1min = pin->GetReal(it->block_name, "x1min");
        ref_size.x1max = pin->GetReal(it->block_name, "x1max");
        if (multi_d) {
          ref_size.x2min = pin->GetReal(it->block_name, "x2min");
          ref_size.x2max = pin->GetReal(it->block_name, "x2max");
        } else {
          ref_size.x2min = mesh_size.x2min;
          ref_size.x2max = mesh_size.x2max;
        }
        if (three_d) {
          ref_size.x3min = pin->GetReal(it->block_name, "x3min");
          ref_size.x3max = pin->GetReal(it->block_name, "x3max");
        } else {
          ref_size.x3min = mesh_size.x3min;
          ref_size.x3max = mesh_size.x3max;
        }
        int phy_ref_lev = pin->GetInteger(it->block_name, "level");
        int log_ref_lev = phy_ref_lev + root_level;
        if (log_ref_lev > current_level) current_level = log_ref_lev;

        // error check parameters in "refinement" blocks
        if (phy_ref_lev < 1) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl <<"<refined_region> level must be larger than 0 (root level=0)"
              << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (log_ref_lev > max_level) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<refined_region> level exceeds maximum allowed ("
              << max_level << ")" << std::endl << "Reduce/specify 'num_levels' in "
              << "<mesh_refinement> input block if using AMR" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (   ref_size.x1min > ref_size.x1max
            || ref_size.x2min > ref_size.x2max
            || ref_size.x3min > ref_size.x3max)  {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Invalid <refined_region> (xmax < xmin in one direction)."
              << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (   ref_size.x1min < mesh_size.x1min || ref_size.x1max > mesh_size.x1max
            || ref_size.x2min < mesh_size.x2min || ref_size.x2max > mesh_size.x2max
            || ref_size.x3min < mesh_size.x3min || ref_size.x3max > mesh_size.x3max) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<refined_region> must be fully contained within root mesh"
              << std::endl;
          std::exit(EXIT_FAILURE);
        }

        // note: if following is too slow, it could be replaced with bi-section search.
        // Suppose entire root domain is tiled with MeshBlocks at the desired refinement
        // level. Find range of x1-integer indices of such MeshBlocks that cover the
        // refinement region
        std::int32_t lx1min = 0, lx1max = 0;
        std::int32_t lx2min = 0, lx2max = 0;
        std::int32_t lx3min = 0, lx3max = 0;
        std::int32_t lxmax = nmb_rootx1*(1<<phy_ref_lev);
        for (lx1min=0; lx1min<lxmax; lx1min++) {
          if (LeftEdgeX(lx1min+1,lxmax,mesh_size.x1min,mesh_size.x1max) > ref_size.x1min)
            break;
        }
        for (lx1max=lx1min; lx1max<lxmax; lx1max++) {
          if (LeftEdgeX(lx1max+1,lxmax,mesh_size.x1min,mesh_size.x1max) >= ref_size.x1max)
            break;
        }
        if (lx1min % 2 == 1) lx1min--;
        if (lx1max % 2 == 0) lx1max++;

        // Find range of x2-indices of such MeshBlocks that cover the refinement region
        if (multi_d) { // 2D or 3D
          lxmax = nmb_rootx2*(1<<phy_ref_lev);
          for (lx2min=0; lx2min<lxmax; lx2min++) {
            if (LeftEdgeX(lx2min+1, lxmax, mesh_size.x2min, mesh_size.x2max) >
                ref_size.x2min)
            break;
          }
          for (lx2max=lx2min; lx2max<lxmax; lx2max++) {
            if (LeftEdgeX(lx2max+1, lxmax, mesh_size.x2min, mesh_size.x2max) >=
                ref_size.x2max)
            break;
          }
          if (lx2min % 2 == 1) lx2min--;
          if (lx2max % 2 == 0) lx2max++;
        }

        // Find range of x3-indices of such MeshBlocks that cover the refinement region
        if (three_d) { // 3D
          lxmax = nmb_rootx3*(1<<phy_ref_lev);
          for (lx3min=0; lx3min<lxmax; lx3min++) {
            if (LeftEdgeX(lx3min+1, lxmax, mesh_size.x3min, mesh_size.x3max) >
                ref_size.x3min)
            break;
          }
          for (lx3max=lx3min; lx3max<lxmax; lx3max++) {
            if (LeftEdgeX(lx3max+1, lxmax, mesh_size.x3min, mesh_size.x3max) >=
                ref_size.x3max)
            break;
          }
          if (lx3min % 2 == 1) lx3min--;
          if (lx3max % 2 == 0) lx3max++;
        }

        // Now add nodes to the MeshBlockTree corresponding to these MeshBlocks
        if (one_d) {  // 1D
          for (std::int32_t i=lx1min; i<lx1max; i+=2) {
            for (int p = 0; p < npanels; ++p) {
              LogicalLocation nlloc;
              nlloc.level = log_ref_lev;
              nlloc.lx1 = i;
              nlloc.lx2 = 0;
              nlloc.lx3 = 0;
              nlloc.panel = p;
              int nnew;
//            ptree->AddNode(nlloc, nnew);
              panel_trees[p]->AddNode(nlloc, nnew);
            }
          }
        }
        if (two_d) {  // 2D
          for (std::int32_t j=lx2min; j<lx2max; j+=2) {
            for (std::int32_t i=lx1min; i<lx1max; i+=2) {
              for (int p = 0; p < npanels; ++p) {
                LogicalLocation nlloc;
                nlloc.level = log_ref_lev;
                nlloc.lx1 = i;
                nlloc.lx2 = j;
                nlloc.lx3 = 0;
                nlloc.panel = p;
                int nnew;
//              ptree->AddNode(nlloc, nnew);
                panel_trees[p]->AddNode(nlloc, nnew);
              }
            }
          }
        }
        if (three_d) {  // 3D
          for (std::int32_t k=lx3min; k<lx3max; k+=2) {
            for (std::int32_t j=lx2min; j<lx2max; j+=2) {
              for (std::int32_t i=lx1min; i<lx1max; i+=2) {
                for (int p = 0; p < npanels; ++p) {
                  LogicalLocation nlloc;
                  nlloc.level = log_ref_lev;
                  nlloc.lx1 = i;
                  nlloc.lx2 = j;
                  nlloc.lx3 = k;
                  nlloc.panel = p;
                  int nnew;
//                ptree->AddNode(nlloc, nnew);
                  panel_trees[p]->AddNode(nlloc, nnew);
                }
              }
            }
          }
        }
      }
    }
  } // if (multilevel)

  if (!adaptive) max_level = current_level;

  // initial mesh hierarchy construction is completed here
//  ptree->CountMeshBlocks(nmb_total);
  nmb_total = 0;
  for (int p = 0; p < npanels; ++p) {
    int nmb_panel;
    panel_trees[p]->CountMeshBlocks(nmb_panel);
    nmb_total += nmb_panel;
  }

  cost_eachmb = new float[nmb_total];
  rank_eachmb = new int[nmb_total];
  lloc_eachmb = new LogicalLocation[nmb_total];
  gids_eachrank = new int[global_variable::nranks];
  nmb_eachrank = new int[global_variable::nranks];

  // following returns LogicalLocation list sorted by Z-ordering, and total # of MBs
//  ptree->CreateZOrderedLLList(lloc_eachmb, nullptr, nmb_total);
  int counttotal = 0;
  for (int p = 0; p < npanels; ++p) {
    panel_trees[p]->CreateZOrderedLLList(counttotal, lloc_eachmb, nullptr, counttotal);
  }

#if MPI_PARALLEL_ENABLED
  // check there is at least one MeshBlock per MPI rank
  if (nmb_total < global_variable::nranks) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
        << "Fewer MeshBlocks (nmb_total=" << nmb_total << ") than MPI ranks (nranks="
        << global_variable::nranks << ")" << std::endl;
    std::exit(EXIT_FAILURE);
  }
#endif

  // initialize cost array with the simplest estimate; all the blocks are equal
  // TODO(@user): implement variable cost per MeshBlock as needed
  for (int i=0; i<nmb_total; i++) {cost_eachmb[i] = 1.0;}
  LoadBalance(cost_eachmb, rank_eachmb, gids_eachrank, nmb_eachrank, nmb_total);

  // create MeshBlockPack for this rank
  int mbp_gids = gids_eachrank[global_variable::my_rank];
  int mbp_gide = mbp_gids + nmb_eachrank[global_variable::my_rank] - 1;
  nmb_thisrank = nmb_eachrank[global_variable::my_rank];

  pmb_pack = new MeshBlockPack(this, mbp_gids, mbp_gide);
  nmb_packs_thisrank = 1;
  pmb_pack->AddMeshBlocks(pin);
//  pmb_pack->pmb->SetNeighbors(ptree, rank_eachmb);
  pmb_pack->pmb->SetNeighbors(rank_eachmb);
  if (use_cubed_sphere && multilevel) {CheckCubedSphereRefinement(pmb_pack);}

  // Fix maximum number of MeshBlocks per rank with AMR
  nmb_maxperrank = nmb_thisrank;
  if (adaptive) {
    if (pin->DoesParameterExist("mesh_refinement", "max_nmb_per_rank")) {
      nmb_maxperrank = pin->GetReal("mesh_refinement", "max_nmb_per_rank");
      if (nmb_maxperrank < nmb_thisrank) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "On rank=" << global_variable::my_rank << " Root grid requires "
          << "more MeshBlocks (nmb_thisrank=" << nmb_thisrank << ") than specified by "
          << "<mesh_refinement>/max_nmb_per_rank=" << nmb_maxperrank << std::endl;
        std::exit(EXIT_FAILURE);
      }
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "With AMR maximum number of MeshBlocks per rank must be "
        << "specified in input file using <mesh_refinement>/max_nmb_per_rank"
        << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
#if MPI_PARALLEL_ENABLED
  if (nmb_maxperrank > (1 << (NUM_BITS_LID))) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
      << "Maximum number of MeshBlocks per rank cannot exceed 2^(NUM_BITS_LID) due to MPI"
      << " tag limits" << std::endl;
    std::exit(EXIT_FAILURE);
  }
#endif

  // Create new MeshRefinement object with either SMR or AMR (SMR needs Restrict fns)
  if (multilevel) {
    pmr = new MeshRefinement(this, pin);
  }

  // set initial time/cycle parameters, output diagnostics
  time = pin->GetOrAddReal("time", "start_time", 0.0);
  dt   = std::numeric_limits<float>::max();
  cfl_no = pin->GetReal("time", "cfl_number");
  ncycle = 0;
  if (global_variable::my_rank == 0) {PrintMeshDiagnostics();}

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::BuildTreeFromRestart():
//! Constructs MeshBlockTree, creates MeshBlockPack (containing the physics modules), and
//! divides grid into MeshBlock(s) for restart runs, using parameters and data read from
//! restart file.

void Mesh::BuildTreeFromRestart(ParameterInput *pin, IOWrapper &resfile,
                                                     bool single_file_per_rank) {
  // At this point, the restartfile is already open and the ParameterInput (input file)
  // data has already been read in main(). Thus the file pointer is set to after <par_end>
  IOWrapperSizeT headeroffset = resfile.GetPosition();

  // following must be identical to calculation of headeroffset (excluding size of
  // ParameterInput data) in restart.cpp
  IOWrapperSizeT headersize = 3*sizeof(int) + 2*sizeof(Real)
    + sizeof(RegionSize) + 2*sizeof(RegionIndcs);
  char *headerdata = new char[headersize];

  // the master process reads the header data if single_file_per_rank is false
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    IOWrapperSizeT read_size = resfile.Read_bytes(headerdata, 1, headersize,
                                                  single_file_per_rank);
    if (read_size != headersize) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Header size read from restart file is incorrect, "
                << "expected " << headersize << ", got " << read_size << std::endl;
      exit(EXIT_FAILURE);
    }
  }

#if MPI_PARALLEL_ENABLED
  // then broadcast the header data
  if (!single_file_per_rank) {
    int mpi_err = MPI_Bcast(headerdata, headersize, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (mpi_err != MPI_SUCCESS) {
      char error_string[1024];
      int length_of_error_string;
      MPI_Error_string(mpi_err, error_string, &length_of_error_string);
      std::cout << "MPI_Bcast failed with error: " << error_string << std::endl;
      exit(EXIT_FAILURE);
    }
  }
#endif

  // Now copy mesh data read from restart file into Mesh variables. Order of variables
  // set by Write()'s in restart.cpp
  // Note this overwrites size and indices initialized in Mesh constructor.
  IOWrapperSizeT hdos = 0;
  std::memcpy(&nmb_total, &(headerdata[hdos]), sizeof(int));
  hdos += sizeof(int);
  std::memcpy(&root_level, &(headerdata[hdos]), sizeof(int));
  hdos += sizeof(int);
  std::memcpy(&mesh_size, &(headerdata[hdos]), sizeof(RegionSize));
  hdos += sizeof(RegionSize);
  std::memcpy(&mesh_indcs, &(headerdata[hdos]), sizeof(RegionIndcs));
  hdos += sizeof(RegionIndcs);
  std::memcpy(&mb_indcs, &(headerdata[hdos]), sizeof(RegionIndcs));
  hdos += sizeof(RegionIndcs);
  std::memcpy(&time, &(headerdata[hdos]), sizeof(Real));
  hdos += sizeof(Real);
  std::memcpy(&dt, &(headerdata[hdos]), sizeof(Real));
  hdos += sizeof(Real);
  std::memcpy(&ncycle, &(headerdata[hdos]), sizeof(int));
  delete [] headerdata;

  // calculate the number of MeshBlocks at root level in each dir
  nmb_rootx1 = mesh_indcs.nx1/mb_indcs.nx1;
  nmb_rootx2 = mesh_indcs.nx2/mb_indcs.nx2;
  nmb_rootx3 = mesh_indcs.nx3/mb_indcs.nx3;
  int current_level = root_level;

  // Error check properties of input paraemters for SMR/AMR meshes.
  if (adaptive) {
    max_level = pin->GetOrAddInteger("mesh_refinement", "num_levels", 1) + root_level - 1;
    if (max_level > 31) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Number of refinement levels must be smaller than "
                << 31 - root_level + 1 << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {
    max_level = 31;
  }

  // allocate memory for lists read from restart
  cost_eachmb = new float[nmb_total];
  rank_eachmb = new int[nmb_total];
  lloc_eachmb = new LogicalLocation[nmb_total];
  gids_eachrank = new int[global_variable::nranks];
  nmb_eachrank = new int[global_variable::nranks];

  // allocate idlist buffer and read list of logical locations and cost
  IOWrapperSizeT listsize = sizeof(LogicalLocation) + sizeof(float);
  char *idlist = new char[listsize*nmb_total];
  // only the master process reads the ID list
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    if (resfile.Read_bytes(idlist,listsize,nmb_total,single_file_per_rank) !=
        static_cast<unsigned int>(nmb_total)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Incorrect number of MeshBlocks in restart file; "
                << "restart file is broken." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
#if MPI_PARALLEL_ENABLED
  // then broadcast the ID list
  if (!single_file_per_rank) {
    MPI_Bcast(idlist, listsize*nmb_total, MPI_CHAR, 0, MPI_COMM_WORLD);
  }
#endif

  // everyone sets the logical location and cost lists based on bradcasted data
  int os = 0;
  for (int i=0; i<nmb_total; i++) {
    std::memcpy(&(lloc_eachmb[i]), &(idlist[os]), sizeof(LogicalLocation));
    os += sizeof(LogicalLocation);
  }
  for (int i=0; i<nmb_total; i++) {
    std::memcpy(&(cost_eachmb[i]), &(idlist[os]), sizeof(float));
    os += sizeof(float);
    if (lloc_eachmb[i].level > current_level) current_level = lloc_eachmb[i].level;
  }
  delete [] idlist;
  if (!adaptive) max_level = current_level;

  // rebuild the MeshBlockTree. Build it through panel_trees, exactly as
  // BuildTreeFromScratch does: that vector owns the nodes and is what the neighbour
  // search in MeshBlock indexes by panel, so a restart that populated only the bare ptree
  // left it empty and every panel lookup reading past the end of it.
  panel_trees.resize(npanels);
  for (int p = 0; p < npanels; ++p) {
    panel_trees[p] = std::make_unique<MeshBlockTree>(this);
    panel_trees[p]->SetPanelID(p);
    panel_trees[p]->CreateRootGrid();
  }
  ptree = panel_trees[0].get();
  // Each block goes into the tree for its OWN panel. The panel id round-trips through the
  // restart file as part of LogicalLocation. Putting them all into panel 0 would be right
  // only for a single-panel mesh: blocks on different panels share logical coordinates,
  // so they would collide.
  for (int i=0; i<nmb_total; i++) {
    int p = lloc_eachmb[i].panel;
    if (p < 0 || p >= npanels) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "MeshBlock " << i << " in the restart file has panel id " << p
        << ", outside the range [0," << npanels-1 << "]. The file was probably written "
        << "by a build with a different number of panels." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    panel_trees[p]->AddNodeWithoutRefinement(lloc_eachmb[i]);
  }

  // check the tree structure by making sure total # of MBs counted in tree same as the
  // number read from the restart file.
  {
    int nnb = 0;
    for (int p = 0; p < npanels; ++p) {
      panel_trees[p]->CreateZOrderedLLList(nnb, lloc_eachmb, nullptr, nnb);
    }
    if (nnb != nmb_total) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "Tree reconstruction failed. Total number of blocks in "
        << "reconstructed tree=" << nnb << ", number in file=" << nmb_total << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

#ifdef MPI_PARALLEL_ENABLED
  // check there is at least one MeshBlock per MPI rank
  if (!single_file_per_rank) {
    if (nmb_total < global_variable::nranks) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line "
        << __LINE__ << std::endl
        << "Fewer MeshBlocks (nmb_total=" << nmb_total << ") than MPI ranks (nranks="
        << global_variable::nranks << ")" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
#endif

  LoadBalance(cost_eachmb, rank_eachmb, gids_eachrank, nmb_eachrank, nmb_total);

  // create MeshBlockPack for this rank
  int mbp_gids = gids_eachrank[global_variable::my_rank];
  int mbp_gide = mbp_gids + nmb_eachrank[global_variable::my_rank] - 1;
  nmb_thisrank = nmb_eachrank[global_variable::my_rank];

  pmb_pack = new MeshBlockPack(this, mbp_gids, mbp_gide);
  pmb_pack->AddMeshBlocks(pin);
//  pmb_pack->pmb->SetNeighbors(ptree, rank_eachmb);
  pmb_pack->pmb->SetNeighbors(rank_eachmb);
  if (use_cubed_sphere && multilevel) {CheckCubedSphereRefinement(pmb_pack);}

  // Fix maximum number of MeshBlocks per rank with AMR
  nmb_maxperrank = nmb_thisrank;
  if (adaptive) {
    if (pin->DoesParameterExist("mesh_refinement", "max_nmb_per_rank")) {
      nmb_maxperrank = pin->GetReal("mesh_refinement", "max_nmb_per_rank");
      if (nmb_maxperrank < nmb_thisrank) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "On rank=" << global_variable::my_rank << " Root grid requires "
          << "more MeshBlocks (nmb_thisrank=" << nmb_thisrank << ") than specified by "
          << "<mesh_refinement>/max_nmb_per_rank=" << nmb_maxperrank << std::endl;
        std::exit(EXIT_FAILURE);
      }
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "With AMR maximum number of MeshBlocks per rank must be "
        << "specified in input file using <mesh_refinement>/max_nmb_per_rank"
        << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // Create new MeshRefinement object with either SMR or AMR (SMR needs Restrict fns)
  if (multilevel) {
    pmr = new MeshRefinement(this, pin);
  }

  // set remaining parameters, output diagnostics
  cfl_no = pin->GetReal("time", "cfl_number");
  if (global_variable::my_rank == 0) {PrintMeshDiagnostics();}
}
