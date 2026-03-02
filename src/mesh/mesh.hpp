#ifndef MESH_MESH_HPP_
#define MESH_MESH_HPP_
//========================================================================================
// Athena++K astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mesh.hpp
//! \brief defines Mesh class.
//! The Mesh is the overall grid structure, which is divided into local patches called
//! MeshBlocks (potentially on different levels) that tile the entire domain.  MeshBlocks
//! are grouped together into MeshBlockPacks for better performance on GPUs.

#include <cstdint>  // int32_t
#include <memory>
#include <string>

#include "athena.hpp"
#include "mesh/nghbr_index.hpp"

// Define following structure before other "include" files to resolve declarations
//----------------------------------------------------------------------------------------
//! \struct RegionSize
//! \brief physical size in a Mesh or a MeshBlock

struct RegionSize {
  Real x1min, x2min, x3min;
  Real x1max, x2max, x3max;
  Real dx1, dx2, dx3;       // (uniform) grid spacing
};

//----------------------------------------------------------------------------------------
//! \struct RegionIndcs
//! \brief Cell indices and number of active and ghost cells in a Mesh or a MeshBlock

struct RegionIndcs {
  int ng;                       // number of ghost cells
  int nx1, nx2, nx3;            // number of active cells (not including ghost zones)
  int is,ie,js,je,ks,ke;        // indices of ACTIVE cells
  int cnx1, cnx2, cnx3;         // number of active coarse cells (not including gzs)
  int cis,cie,cjs,cje,cks,cke;  // indices of ACTIVE coarse cells
};

//----------------------------------------------------------------------------------------
//! \struct NeighborBlock
//! \brief Information about neighboring MeshBlocks stored as 2D DualArray in MeshBlock

struct NeighborBlock {
  int gid;     // global ID
  int lev;     // logical level
  int rank;    // MPI rank
  int dest;    // index of recv buffer in target NeighborBlocks
  int panel;   // panel ID
};

//----------------------------------------------------------------------------------------
//! \struct LogicalLocation
//! \brief logical location and level of MeshBlock stored as POD
//! lx1/2/3 = logical location in x1/2/3 = index in array of nodes at current level
//! WARNING: values of lx? can exceed the range of std::int32_t with >30 levels
//! of AMR, even if the root grid consists of a single MeshBlock, since the corresponding
//! max index = 1*2^31 > INT_MAX = 2^31 -1 for most 32-bit signed integer types

struct LogicalLocation {
  std::int32_t lx1, lx2, lx3, level, panel;
};

//----------------------------------------------------------------------------------------
//! \struct EventCounters
//! \brief stores various counters used as diagnostics throughout the code

struct EventCounters {
  int nfofc, neos_dfloor, neos_efloor, neos_tfloor, neos_vceil, neos_fail, maxit_c2p;
  EventCounters() : nfofc(0), neos_dfloor(0), neos_efloor(0), neos_tfloor(0),
                    neos_vceil(0), neos_fail(0), maxit_c2p(0) {}
};

//----------------------------------------------------------------------------------------
//! \struct PanelNeighbors
//! \brief stores connectivity table
//!
struct PanelNeighbors {
  int nb_panel;   // neighbor panel id (0-6)
  int nb_face;    // which face of neighbor we enter (0:-x1, 1:+x1, 2:-x2, 3:+x2)
  int rev_ax;     // whether along-edge index is reversed (0/1)
  int swap_ax;    // whether to swap the xy axis (0/1)
};

struct PanelBoundaries {
  int swap_ax;    // whether to swap the xy axis (0/1)
  int rev_x1;     // whether (after being swapped) x1 index is reversed (0/1)
  int rev_x2;     // whether (after being swapped) x2 index is reversed (0/1)
  int end_face;   // which face of the panel we transform into (0:-x1, 1:+x1, 2:-x2, 3:+x2)
};

// Forward declarations required due to recursive definitions amongst mesh classes
class MeshBlock;
class MeshBlockPack;
class MeshBlockTree;
class Mesh;

#include "parameter_input.hpp"
#include "meshblock.hpp"
#include "meshblock_pack.hpp"
#include "meshblock_tree.hpp"
#include "mesh_refinement.hpp"

//----------------------------------------------------------------------------------------
//! \class Mesh
//! \brief data/functions associated with the overall mesh

class Mesh {
  // mesh classes (Mesh, MeshBlock, MeshBlockPack, MeshBlockTree, MeshRefinement)
  // like to play together
  friend class MeshBlock;
  friend class MeshBlockPack;
  friend class MeshBlockTree;
  friend class MeshRefinement;
  // needs to access tree to find target MB offset by shear
  friend class ShearingBox;

 public:
  explicit Mesh(ParameterInput *pin);
  ~Mesh();

  // data
  RegionSize  mesh_size;      // (physical) size of mesh (physical root level)
  RegionIndcs mesh_indcs;     // indices of cells in mesh (physical root level)
  RegionIndcs mb_indcs;       // indices of cells in MeshBlocks (same for all MeshBlocks)
  BoundaryFlag mesh_bcs[6];   // physical boundary conditions at 6 faces of mesh
  bool strictly_periodic;     // true if all boundaries are periodic
    
  bool use_cubed_sphere;      // true if using cubed sphere
  int npanels;                // 6 if using cubed sphere; 1 otherwise

  bool one_d, two_d, three_d; // flags to indicate 1D or 2D or 3D calculations
  bool multi_d;               // flag to indicate 2D and 3D calculations
  bool multilevel;            // true for SMR and AMR
  bool adaptive;              // true only for AMR

  int nmb_rootx1, nmb_rootx2, nmb_rootx3; // # of MeshBlocks at root level in each dir
  int nmb_total;           // total number of MeshBlocks across all levels/ranks
  int nmb_thisrank;        // number of MeshBlocks on this MPI rank (local)
  int nmb_maxperrank;      // max allowed number of MBs per device (memory limit for AMR)

  int root_level; // logical level of root (physical) grid (e.g. Fig. 3 of method paper)
  int max_level;  // logical level of maximum refinement grid in Mesh

  int nprtcl_thisrank;     // number of particles this rank
  int nprtcl_total;        // total number of particles across all ranks

  // following 3x arrays allocated with length [nmb_total] in BuildTreeFromXXXX()
  float *cost_eachmb;            // cost of each MeshBlock
  int *rank_eachmb;              // rank of each MeshBlock
  LogicalLocation *lloc_eachmb;  // LogicalLocations for each MeshBlock

  // following 2x arrays allocated with length [nranks] in BuildTreeFromXXXX()
  int *gids_eachrank;      // starting global ID of MeshBlocks in each rank
  int *nmb_eachrank;       // number of MeshBlocks on each rank
  // following 1x arrays allocated with length [nranks] in AddCoordinatesAndPhysics()
  int *nprtcl_eachrank;    // number of particles on each rank

  Real time, dt, dtold, cfl_no;
  int ncycle;
  EventCounters ecounter;

  int nmb_packs_thisrank;                  // number of MBPacks on this rank
  MeshBlockPack* pmb_pack;                 // container for MeshBlocks on this rank
  std::unique_ptr<ProblemGenerator> pgen;  // class containing functions to set ICs
  MeshRefinement *pmr=nullptr;             // mesh refinement data/functions (if needed)

  // functions
  void BuildTreeFromScratch(ParameterInput *pin);
  void BuildTreeFromRestart(ParameterInput *pin, IOWrapper &resfile,
                            bool single_file_per_rank=false);
  void PrintMeshDiagnostics();
  void WriteMeshStructure();
  void NewTimeStep(const Real tlim);
  void AddCoordinatesAndPhysics(ParameterInput *pinput);
  BoundaryFlag GetBoundaryFlag(const std::string& input_string);
  std::string GetBoundaryString(BoundaryFlag input_flag);
    
    KOKKOS_INLINE_FUNCTION
    int NeighborIndexPanel(int ix, int iy, int iz, int n1, int n2, int panel_start, int panel_end) {
        
        const PanelBoundaries pb = GetPanelBoundary(panel_start, panel_end);
        int iix = (pb.swap_ax == 1) ? iy : ix;
        int iiy = (pb.swap_ax == 1) ? ix : iy;
        if (pb.rev_x1 == 1) iix = -iix;
        if (pb.rev_x2 == 1) iiy = -iiy;
        
        return NeighborIndex(iix,iiy,iz,n1,n2);
    }

  // comparison function for sorting LogicalLocations based on level
  static bool GreaterLevel(const LogicalLocation & left, const LogicalLocation &right) {
    return left.level > right.level;
  }

  // accessors
  int FindMeshBlockIndex(int tgid) {
    for (int m=0; m<pmb_pack->nmb_thispack; ++m) {
      if (pmb_pack->pmb->mb_gid.h_view(m) == tgid) return m;
    }
    return -1;
  }
  int NumberOfMeshBlockCells() const {
    return (mb_indcs.nx1)*(mb_indcs.nx2)*(mb_indcs.nx3);
  }
    
    /* Cubed sphere. Borrowed from SNAPY https://github.com/chengcli/snapy/blob/main/src/layout/cubed_sphere_layout.cpp
    *                                            z  y
    *                 ___________                | /
    *                 |\        .\               |/---x
    *                 | \   3   . \             (3)
    *                 |  \_________\
    *      y          | 4 |     .  |
    *      |          \. .|......  |
    *  z___|(4)        \  |    0 . |
    *     /             \ |       .|             y
    *   x/               \|________|             |
    *                                        (0) |___x
    *                                           /
    *                                         z/
    *
    *
    *
    *                                         y  x
    *                 __________              | /
    *                 |\       .\             |/___z
    *                 | \      . \           (1)
    *      y  z       |  \________\
    *      | /        |  |  2  .  |
    *  x___|/         |..|......  |
    *       (2)       \  |     . 1|        (5)  ___ x
    *                  \ |  5   . |           /|
    *                   \|_______.|         y/ |
    *                                          z
    *
    * --------------------------
    * Cubed-sphere connectivity
    * --------------------------
    * Face numbering (editable):
    *
    *           -------
    *           |  3  |
    *     |-----|-----|-----|-----|
    *     |  4  |  0  |  1  |  2  |
    *     |-----|-----|-----|-----|
    *           |  5  |
    *           |-----|
    *
    */
    
    // Face ordering:
    // 0:-x1, 1:+x1, 2:-x2, 3:+x2, 4:-x3, 5:+x3

//    static constexpr PanelNeighbors panel_neighbors[6][4] = {
//      // panel 0: neighbors 4(L),1(R),5(B),3(T)
//      {
//        {4, 1, 0, 0},  // L
//        {1, 0, 0, 0},  // R
//        {5, 3, 0, 0},  // B
//        {3, 2, 0, 0}   // T
//      },
//
//      // panel 1: neighbors 0(L),2(R),5(B),3(T)
//      {
//        {0, 1, 0, 0},  // L
//        {2, 0, 0, 0},  // R
//        {5, 1, 1, 1},  // B
//        {3, 1, 0, 1}   // T
//      },
//
//      // panel 2: neighbors 1(L),4(R),5(B),3(T)
//      {
//        {1, 1, 0, 0},  // L
//        {4, 0, 0, 0},  // R
//        {5, 2, 1, 0},  // B
//        {3, 3, 1, 0}   // T
//      },
//
//      // panel 3: neighbors 4(L),1(R),0(B),2(T)
//      {
//        {4, 3, 1, 1},  // L
//        {1, 3, 0, 1},  // R
//        {0, 3, 0, 0},  // B
//        {2, 3, 1, 0}   // T
//      },
//
//      // panel 4: neighbors 2(L),0(R),5(B),3(T)
//      {
//        {2, 1, 0, 0},  // L
//        {0, 0, 0, 0},  // R
//        {5, 0, 0, 1},  // B
//        {3, 0, 1, 1}   // T
//      },
//
//      // panel 5: neighbors 4(L),1(R),2(B),0(T)
//      {
//        {4, 2, 0, 1},  // L
//        {1, 2, 1, 1},  // R
//        {2, 2, 1, 0},  // B
//        {0, 2, 0, 0}   // T
//      }
//    };
    
    static constexpr PanelNeighbors panel_neighbors[2][4] = {
      // panel 0
      {
        {0, 1, 0, 0},  // L
        {0, 0, 0, 0},  // R
        {1, 0, 0, 1},  // B
        {1, 1, 0, 1},  // T
      },

      // panel 1
      {
        {0, 2, 0, 1}, // L
        {0, 3, 0, 1}, // R
        {1, 3, 0, 0}, // B
        {1, 2, 0, 0}  // T
      }
    };
    
    KOKKOS_INLINE_FUNCTION
    constexpr PanelBoundaries GetPanelBoundary(int ps, int pe) const {
      constexpr PanelBoundaries table[2][2] = {
        // … your SAME initializer …
          // [i][j]: transforming from i panel to j panel
          // panel 0:
          {
            {0, 0, 0, -1}, // 00
            {1, 1, 0, 0},  // 01
          },

          // panel 1:
          {
            {1, 0, 1, 0},  // 10
            {0, 0, 0, -1}, // 11
          }
      };
      return table[ps][pe];
    }
    
//    KOKKOS_INLINE_FUNCTION
//    constexpr PanelBoundaries GetPanelBoundary(int ps, int pe) const {
//      constexpr PanelBoundaries table[6][6] = {
//        // … your SAME initializer …
//          // [i][j]: transforming from i panel to j panel
//          // panel 0:
//          {
//            {0, 0, 0, -1}, // 00
//            {0, 0, 0, 0},  // 01
//            {0, 0, 0, -1}, // 02
//            {0, 0, 0, 2},  // 03
//            {0, 0, 0, 1},  // 04
//            {0, 0, 0, 3}   // 05
//          },
//
//          // panel 1:
//          {
//            {0, 0, 0, 1},  // 10
//            {0, 0, 0, -1}, // 11
//            {0, 0, 0, 0},  // 12
//            {1, 1, 0, 1},  // 13
//            {0, 0, 0, -1}, // 14
//            {1, 0, 1, 1}   // 15
//          },
//
//          // panel 2:
//          {
//            {0, 0, 0, -1}, // 20
//            {0, 0, 0, 1},  // 21
//            {0, 0, 0, -1}, // 22
//            {0, 1, 1, 3},  // 23
//            {0, 0, 0, 0},  // 24
//            {0, 1, 1, 2}   // 25
//          },
//
//          // panel 3:
//          {
//            {0, 0, 0, 3},  // 30
//            {1, 0, 1, 3},  // 31
//            {0, 1, 1, 3},  // 32
//            {0, 0, 0, -1}, // 33
//            {1, 1, 0, 3},  // 34
//            {0, 0, 0, -1}  // 35
//          },
//
//          // panel 4:
//          {
//            {0, 0, 0, 0},  // 40
//            {0, 0, 0, -1}, // 41
//            {0, 0, 0, 1},  // 42
//            {1, 0, 1, 0},  // 43
//            {0, 0, 0, -1}, // 44
//            {1, 1, 0, 0}   // 45
//          },
//
//          // panel 5:
//          {
//            {0, 0, 0, 2},  // 50
//            {1, 1, 0, 2},  // 51
//            {0, 1, 1, 2},  // 52
//            {0, 0, 0, -1}, // 53
//            {1, 0, 1, 2},  // 54
//            {0, 0, 0, -1}   // 55
//          }
//      };
//      return table[ps][pe];
//    }


 private:
  std::unique_ptr<MeshBlockTree> ptree;  // pointer to root node in binary/quad/oct-tree
  std::vector<std::unique_ptr<MeshBlockTree>> panel_trees;
  void LoadBalance(float *clist, int *rlist, int *slist, int *nlist, int nb);
};
#endif  // MESH_MESH_HPP_
