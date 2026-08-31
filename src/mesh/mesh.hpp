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
#include "coordinates/grid_stretch.hpp"
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
  int polar;   // 1 = polar boundary; -1 = otherwise
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
// NOTE ON AXES. On the cubed sphere x1 is the RADIAL direction, matching
// CoordSphericalPolar, and the two PANEL-TANGENTIAL axes are x2 (xi) and x3 (eta). A
// panel seam therefore permutes and reflects x2 and x3 only; x1 is never touched by
// one. The two tangential axes are named `a` (= x2) and `b` (= x3) throughout the
// panel tables below, deliberately NOT `x1`/`x2`: they used to be, under an earlier
// convention in which the radius was x3, and reading those names as mesh axes is what
// silently reversed the radial index across every seam.
// Face ordering everywhere in the panel tables: 0:-x2, 1:+x2, 2:-x3, 3:+x3.

struct PanelNeighbors {
  int nb_panel;   // neighbor panel id (0-6)
  int nb_face;    // which face of neighbor we enter (0:-x2, 1:+x2, 2:-x3, 3:+x3)
  int rev_ax;     // whether along-edge index is reversed (0/1)
  int swap_ax;    // whether to swap the two tangential axes a,b (0/1)
};

struct PanelBoundaries {
  int swap_ax;    // whether to swap the two tangential axes a,b (0/1)
  int rev_a;      // whether (after being swapped) the a = x2 index is reversed (0/1)
  int rev_b;      // whether (after being swapped) the b = x3 index is reversed (0/1)
  int end_face;   // which face of the starting panel (0:-x2, 1:+x2, 2:-x3, 3:+x3)
};

//----------------------------------------------------------------------------------------
//! \fn GetPanelBoundary
//! \brief How panel `ps` meets panel `pe`: whether the two tangential axes swap, and
//! whether either index reverses.
//!
//! A FREE function, not a member of Mesh, and deliberately so: every caller is a device
//! kernel, and reaching it through `pmy_pack->pmesh->` made the lambda capture `this` and
//! dereference a host pointer on the device -- which hipcc warns about
//! (-Wgpu-maybe-wrong-side) and which is an illegal access on any GPU that does not share
//! memory with the host. It reads nothing but its two arguments, so it never needed to be
//! a member. See also cubed_sphere::PanelTangents, which is a free function for the same
//! reason.
KOKKOS_INLINE_FUNCTION
constexpr PanelBoundaries GetPanelBoundary(int ps, int pe) {
  constexpr PanelBoundaries table[6][6] = {
    // … your SAME initializer …
      // [i][j]: transforming from i panel to j panel
      // panel 0:
      {
        {0, 0, 0, -1}, // 00
        {0, 0, 0, 1},  // 01
        {0, 0, 0, -1}, // 02
        {0, 0, 0, 3},  // 03
        {0, 0, 0, 0},  // 04
        {0, 0, 0, 2}   // 05
      },

      // panel 1:
      {
        {0, 0, 0, 0},  // 10
        {0, 0, 0, -1}, // 11
        {0, 0, 0, 1},  // 12
        {1, 1, 0, 3},  // 13
        {0, 0, 0, -1}, // 14
        {1, 0, 1, 2}   // 15
      },

      // panel 2:
      {
        {0, 0, 0, -1}, // 20
        {0, 0, 0, 0},  // 21
        {0, 0, 0, -1}, // 22
        {0, 1, 1, 3},  // 23
        {0, 0, 0, 1},  // 24
        {0, 1, 1, 2}   // 25
      },

      // panel 3:
      {
        {0, 0, 0, 2},  // 30
        {1, 0, 1, 1},  // 31
        {0, 1, 1, 3},  // 32
        {0, 0, 0, -1}, // 33
        {1, 1, 0, 0},  // 34
        {0, 0, 0, -1}  // 35
      },

      // panel 4:
      {
        {0, 0, 0, 1},  // 40
        {0, 0, 0, -1}, // 41
        {0, 0, 0, 0},  // 42
        {1, 0, 1, 3},  // 43
        {0, 0, 0, -1}, // 44
        {1, 1, 0, 2}   // 45
      },

      // panel 5:
      {
        {0, 0, 0, 3},  // 50
        {1, 1, 0, 1},  // 51
        {0, 1, 1, 2},  // 52
        {0, 0, 0, -1}, // 53
        {1, 0, 1, 0},  // 54
        {0, 0, 0, -1}   // 55
      }
  };
  return table[ps][pe];
}

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
//! \struct PanelEdgeMap
//! \brief How a cubed-sphere EDGE direction looks from the neighbouring panel.
//!
//! Everything a mixed-level cross-panel diagonal needs and did not have. Three things
//! differ across a seam and all three were previously ignored at a level boundary:
//!
//!   * the OFFSETS pointing back at this block, which decide which of the neighbour's
//!     children touch it (`GetLeaf` is indexed in the NEIGHBOUR's axes, not ours);
//!   * WHICH of the neighbour's axes is the free one along the edge -- a seam can send
//!     our x2 to their x3, so enumerating their children along the wrong axis picks the
//!     wrong two blocks entirely. That is what made a coarse block list the two fine
//!     blocks stacked in x3 when its actual neighbours were the two stacked in x2, and
//!     it is the reason those slots did not pair;
//!   * whether that free axis is REVERSED, which swaps the two subblock halves.
//!
//! Reduces to the identity when the two panels share a chart, so the same-panel path is
//! untouched.

struct PanelEdgeMap {
  int ox1r, ox2r, ox3r;   // offsets in the NEIGHBOUR's axes, pointing back at us
  int free_axis;          // which of the NEIGHBOUR's axes runs ALONG the edge (1, 2, 3)
  int rev_free;           // 1 if our free axis maps to a reversed neighbour axis
};

KOKKOS_INLINE_FUNCTION
PanelEdgeMap GetPanelEdgeMap(int ox1, int ox2, int ox3, int ps, int pe) {
  const PanelBoundaries pb = GetPanelBoundary(ps, pe);
  PanelEdgeMap pm;
  int ia = (pb.swap_ax == 1) ? (-ox3) : (-ox2);
  int ib = (pb.swap_ax == 1) ? (-ox2) : (-ox3);
  if (pb.rev_a == 1) ia = -ia;
  if (pb.rev_b == 1) ib = -ib;
  pm.ox1r = -ox1;   pm.ox2r = ia;   pm.ox3r = ib;
  pm.free_axis = 1;                          // x2x3 edge: the edge runs along x1
  pm.rev_free = 0;                           // and x1 is radial, which no seam crosses
  if (ox2 == 0 && ox3 != 0) {                // x3x1 edge: our free axis is x2
    pm.free_axis = (pb.swap_ax == 1) ? 3 : 2;
    pm.rev_free = (pb.swap_ax == 1) ? pb.rev_b : pb.rev_a;
  } else if (ox3 == 0 && ox2 != 0) {         // x1x2 edge: our free axis is x3
    pm.free_axis = (pb.swap_ax == 1) ? 2 : 3;
    pm.rev_free = (pb.swap_ax == 1) ? pb.rev_a : pb.rev_b;
  }
  return pm;
}

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
  bool cs_vertex_fill;        // true to fill the cube-vertex corner by exchange
  int npanels;                // 6 if using cubed sphere; 1 otherwise
  bool use_spherical_polar;   // true if using spherical polar grid
  bool use_grid_stretch_r;      // true if using grid stretching in r
  bool use_grid_stretch_r_poly; // true if using the polynomial radial stretch
  bool use_grid_stretch_theta;  // true if using grid stretching in theta
  Real fStretchR, fStretchTheta;
  // Coefficients of the polynomial radial stretch. See StretchRPoly in coordinates.hpp.
  Real fStretchRPoly[NSTRETCH_R_POLY];
  bool use_polar_boundary;      // true if using polar boundaries

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
  Real dt_diff;
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
        
        // ix is the RADIAL (x1) offset and is passed through untouched: a panel seam
        // never crosses it. Only the tangential pair a = x2 (iy), b = x3 (iz) is
        // remapped.
        const PanelBoundaries pb = GetPanelBoundary(panel_start, panel_end);
        int iia = (pb.swap_ax == 1) ? iz : iy;
        int iib = (pb.swap_ax == 1) ? iy : iz;
        if (pb.rev_a == 1) iia = -iia;
        if (pb.rev_b == 1) iib = -iib;
        
        return NeighborIndex(ix,iia,iib,n1,n2);
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
    
    // Face ordering (panel-tangential faces only; a = x2, b = x3):
    // 0:-x2, 1:+x2, 2:-x3, 3:+x3

    static constexpr PanelNeighbors panel_neighbors[6][4] = {
      // panel 0: neighbors 4(L),1(R),5(B),3(T)
      {
        {4, 1, 0, 0},  // L
        {1, 0, 0, 0},  // R
        {5, 3, 0, 0},  // B
        {3, 2, 0, 0}   // T
      },

      // panel 1: neighbors 0(L),2(R),5(B),3(T)
      {
        {0, 1, 0, 0},  // L
        {2, 0, 0, 0},  // R
        {5, 1, 1, 1},  // B
        {3, 1, 0, 1}   // T
      },

      // panel 2: neighbors 1(L),4(R),5(B),3(T)
      {
        {1, 1, 0, 0},  // L
        {4, 0, 0, 0},  // R
        {5, 2, 1, 0},  // B
        {3, 3, 1, 0}   // T
      },

      // panel 3: neighbors 4(L),1(R),0(B),2(T)
      {
        {4, 3, 1, 1},  // L
        {1, 3, 0, 1},  // R
        {0, 3, 0, 0},  // B
        {2, 3, 1, 0}   // T
      },

      // panel 4: neighbors 2(L),0(R),5(B),3(T)
      {
        {2, 1, 0, 0},  // L
        {0, 0, 0, 0},  // R
        {5, 0, 0, 1},  // B
        {3, 0, 1, 1}   // T
      },

      // panel 5: neighbors 4(L),1(R),2(B),0(T)
      {
        {4, 2, 0, 1},  // L
        {1, 2, 1, 1},  // R
        {2, 2, 1, 0},  // B
        {0, 2, 0, 0}   // T
      }
    };
    
//    static constexpr PanelNeighbors panel_neighbors[2][4] = {
//      // panel 0
//      {
//        {0, 1, 0, 0},  // L
//        {0, 0, 0, 0},  // R
//        {1, 0, 0, 1},  // B
//        {1, 1, 0, 1},  // T
//      },
//
//      // panel 1
//      {
//        {0, 2, 0, 1}, // L
//        {0, 3, 0, 1}, // R
//        {1, 3, 0, 0}, // B
//        {1, 2, 0, 0}  // T
//      }
//    };
//
//    KOKKOS_INLINE_FUNCTION
//    constexpr PanelBoundaries GetPanelBoundary(int ps, int pe) const {
//      constexpr PanelBoundaries table[2][2] = {
//        // … your SAME initializer …
//          // [i][j]: transforming from i panel to j panel
//          // panel 0:
//          {
//            {0, 0, 0, -1}, // 00
//            {1, 1, 0, 0},  // 01
//          },
//
//          // panel 1:
//          {
//            {1, 0, 1, 0},  // 10
//            {0, 0, 0, -1}, // 11
//          }
//      };
//      return table[ps][pe];
//    }
    

 private:
  // The panel trees OWN every node. There is one panel for an ordinary mesh and six for a
  // cubed sphere, each with its own root and its own z-ordering.
  std::vector<std::unique_ptr<MeshBlockTree>> panel_trees;
  // Non-owning pointer to the root of the first panel tree. Most of the code predates the
  // panel decomposition and reaches the tree through this; it is the whole tree unless
  // the mesh is a cubed sphere. It MUST be set wherever panel_trees is built -- leaving
  // it null is what made adaptive refinement segfault in MeshBlockTree::FindMeshBlock.
  MeshBlockTree *ptree = nullptr;
  void LoadBalance(float *clist, int *rlist, int *slist, int *nlist, int nb);
};
#endif  // MESH_MESH_HPP_
