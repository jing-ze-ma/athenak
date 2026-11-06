///======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in the code
 * distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

// C++ headers
#include <iostream>   // endl
#include <fstream>    
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>      // sqrt, cos
#include <algorithm>  // min
#include <cstdio>     // printf()
#include <cstdlib>    // srand

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../bvals/bvals.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"
#include "../nr_radiation/radiation.hpp"
#include "../nr_radiation/implicit/radiation_implicit.hpp"
#include "../nr_radiation/integrators/rad_integrators.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "../utils/utils.hpp"
#include "../globals.hpp"
#include "../units/units.hpp" // to get grav_const_cgs, k_boltzmann_cgs, hydrogen_mass_cgs

// Units, need to set manually

static Real rho_unit;
static Real T_unit;
static Real l_unit;
static Real time_unit;
static Real grav_const_dimensionless;

// The space for opacity table, for tables, rows are different temperature and columns are different density

static AthenaArray<Real> log_kappa_ross_table;
static AthenaArray<Real> log_kappa_planck_table;
static AthenaArray<Real> log_kappa_electron_table;
static AthenaArray<Real> log_T_table;
static AthenaArray<Real> log_rho_table;

static int total_num_log_rho_grid;
static int total_num_log_T_grid;

// for initial profile and bdr conditions

static AthenaArray<Real> ini_profile;   // the initial MESA profile: 0: radius, 1: temperature, 2: density, 3: mass, 4: kappa, 5: Fr
static AthenaArray<Real> bd_data;       // the array to save the bdr values, to be assigned to ghost cells
static int num_of_points_init;          // num of points in the initial MESA profile


// some global running variables
static Real tfloor;       // temperature floor value
static Real rhofloor;     // density floor value

static int n_user_var=4;  // num of user output

// This is some function/variables needed to construct face states to pass to boundary conditions
// so that the net mass flux can be forced to zero

struct X1FaceStates {
  Real rho_l, vel_l, pressure_l, sound_speed_l;
  Real rho_r, vel_r, pressure_r, sound_speed_r;
};

void left_and_right_states_at_boundary_hllc(
      MeshBlock *pmb, const AthenaArray<Real> &prim,
      const AthenaArray<Real> &bcc, int k, int j, int is, int ie,
      AthenaArray<Real> &wl, AthenaArray<Real> &wr, X1FaceStates &states);

//======================================================================================
/*! \file bsg_star.cpp
 *  \brief blue supergiant stellar envelope with radiative inner bdr (fixed flux, temperature, density, zero radial velocity)
 *
 *====================================================================================*/

void Inflow_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void Inflow_rad_X1(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
     const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void Outflow_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void Outflow_rad_X1(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
     const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void StarOpacity(MeshBlock *pmb, AthenaArray<Real> &prim);

void interpolate_opacity(const Real rho_cgs, const Real T_gas_cgs, Real &kappa_ross, Real &kappa_planck, Real &kappa_electron);



// col4x4: his 256 radial faces as explicit piecewise-uniform spacing (no SMR)
// 48 x 0.3125 (20-35), 128 x 0.15625 (35-55), 80 x 0.3125 (55-80); linear past the ends
static Real HisGridX1(Real x, RegionSize rs) {
  Real n = x*256.0;
  if (n <= 48.0) return 20.0 + n*0.3125;
  if (n <= 176.0) return 35.0 + (n - 48.0)*0.15625;
  return 55.0 + (n - 176.0)*0.3125;
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  if (mesh_size.x1rat < 0.0) EnrollUserMeshGenerator(X1DIR, HisGridX1);
  
  // enroll boundary functions
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, Outflow_X1);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, Inflow_X1);

  // get hydro and radiation floor values
  tfloor = pin->GetReal("radiation", "tfloor");
  rhofloor = pin->GetReal("hydro", "dfloor");

  // get units
  rho_unit = pin->GetReal("radiation", "density_unit");
  T_unit = pin->GetReal("radiation", "T_unit");
  l_unit = pin->GetReal("radiation", "length_unit");
  Real mu_molecular = pin->GetReal("radiation", "molecular_weight");

  Real sound_speed_cgs = sqrt(Constants::k_boltzmann_cgs * T_unit / (mu_molecular * Constants::hydrogen_mass_cgs));
  time_unit = l_unit / sound_speed_cgs;

  grav_const_dimensionless = Constants:: grav_const_cgs / (1.0/(rho_unit * time_unit * time_unit));

  // enroll gravity potential function, which is the a source term
  // EnrollUserExplicitSourceFunction(GravityPotential);

  // open initial MESA profile
  FILE *fini;
  if ( (fini=fopen("./MESA_profile_combined.txt","r"))==NULL )
  {
     printf("Open input file error MESA profile");
     return;
  }

  // count how many lines in the initial profile

  FILE *fini_countline; // we don't want to reach end of the file for fini, so we reread it solely for line counting
  num_of_points_init = 0;

  fini_countline = fopen("./MESA_profile_combined.txt","r");
  char buffer[2048]; // an empty line of string to pass to fgets function

  while (fgets(buffer, sizeof(buffer), fini_countline)) {
      ++num_of_points_init;
  }

  // read initial radial profile from file, 0: radius, 1: temperature, 2: density, 3: mass, 4: kappa, 5: Fr
  ini_profile.NewAthenaArray(num_of_points_init,6);

  for(int j=0; j<num_of_points_init; j++){
    for(int i=0; i<6; i++){
      fscanf(fini,"%lf",&(ini_profile(j,i)));
    }
  }

  fclose(fini);
  
  // get maximum radius in the initial profile
  Real ini_profile_max_r = ini_profile(num_of_points_init-1,0);

  // inner boundary data array: density, temperature, Fr, in each ghost cell
  bd_data.NewAthenaArray(3,NGHOST);


  // total number of face points, attach NGHOST ghost cells in both directions to account for the bdr conditions
  //---------------------------- e.g. nx1 = 11, NGHOST = 2 -------------------------//
  //|GHOST|GHOST|CELL|CELL|CELL|CELL|CELL|CELL|CELL|CELL|CELL|CELL|CELL|GHOST|GHOST|//
  //0     1     2    3    4    5    6    7    8    9    10   11   12   13    14   15//
  //          x1min                                                                 //
  //--------------------------------------------------------------------------------//


  // get x1min from input, see explanations there
  Real x1min = mesh_size.x1min;

  if(NR_RADIATION_ENABLED || IM_RADIATION_ENABLED){

    // enroll radiative transfer bdr conditions
  
    EnrollUserRadBoundaryFunction(BoundaryFace::inner_x1, Inflow_rad_X1);
    EnrollUserRadBoundaryFunction(BoundaryFace::outer_x1, Outflow_rad_X1);
    
    // read the opacity table files
    FILE *f_log_kappa_ross, *f_log_kappa_planck, *f_log_kappa_electron, *f_log_T, *f_log_rho;

    if ( (f_log_kappa_ross = fopen("./log_kappa_ross.txt","r")) == NULL )
    {
      printf("cannot find log_kappa_ross.txt");
      return;
    }

    if ( (f_log_kappa_planck = fopen("./log_kappa_planck.txt","r")) == NULL )
    {
      printf("cannot find log_kappa_planck.txt");
      return;
    }

    if ( (f_log_kappa_electron = fopen("./log_kappa_electron.txt","r")) == NULL )
    {
      printf("cannot find log_kappa_electron.txt");
      return;
    }

    if ( (f_log_T = fopen("./log_T.txt","r")) == NULL )
    {
      printf("cannot find log_T.txt");
      return;
    }

    if ( (f_log_rho = fopen("./log_rho.txt","r")) == NULL )
    {
      printf("cannot find log_rho.txt");
      return;
    }

    // count the dimensions of log_T and log_rho grids

    FILE *f_log_T_countline, *f_log_rho_countline; // we don't want to reach end of the file, so we reread it solely for line counting
    total_num_log_rho_grid = 0;
    total_num_log_T_grid = 0;

    f_log_T_countline = fopen("./log_T.txt","r");
    f_log_rho_countline = fopen("./log_rho.txt","r");

    while (fgets(buffer, sizeof(buffer), f_log_T_countline)) {
        ++total_num_log_T_grid;
    }

    while (fgets(buffer, sizeof(buffer), f_log_rho_countline)) {
        ++total_num_log_rho_grid;
    }

    // initialize arrays

    log_kappa_ross_table.NewAthenaArray(total_num_log_T_grid,total_num_log_rho_grid);
    log_kappa_planck_table.NewAthenaArray(total_num_log_T_grid,total_num_log_rho_grid);
    log_kappa_electron_table.NewAthenaArray(total_num_log_T_grid,total_num_log_rho_grid);
    log_T_table.NewAthenaArray(total_num_log_T_grid);
    log_rho_table.NewAthenaArray(total_num_log_rho_grid);

    // assign table values to Athena arrays

    for(int j = 0; j < total_num_log_T_grid; j++){
      for(int i = 0; i < total_num_log_rho_grid; i++){
          fscanf(f_log_kappa_ross, "%lf", &(log_kappa_ross_table(j,i)));
          fscanf(f_log_kappa_planck, "%lf", &(log_kappa_planck_table(j,i)));
          fscanf(f_log_kappa_electron, "%lf", &(log_kappa_electron_table(j,i)));
      }
    }

    for(int i = 0; i < total_num_log_T_grid; i++){
      fscanf(f_log_T,"%lf",&(log_T_table(i)));
    }

    for(int i = 0; i < total_num_log_rho_grid; i++){
      fscanf(f_log_rho,"%lf",&(log_rho_table(i)));
    }

    // close files
    fclose(f_log_kappa_ross);
    fclose(f_log_kappa_electron);
    fclose(f_log_kappa_planck);
    fclose(f_log_rho);
    fclose(f_log_T);

  }

  return;
}

void Mesh::UserWorkInLoop()
{
  // nothing needs to be done here
}
//======================================================================================
//! \fn void Mesh::TerminateUserMeshProperties(void)
//  \brief Clean up the Mesh properties
//======================================================================================
void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{

  ini_profile.DeleteAthenaArray();
  bd_data.DeleteAthenaArray();

  // free memory
  if(NR_RADIATION_ENABLED || IM_RADIATION_ENABLED){

    log_kappa_ross_table.DeleteAthenaArray();
    log_kappa_planck_table.DeleteAthenaArray();
    log_kappa_electron_table.DeleteAthenaArray();
    log_T_table.DeleteAthenaArray();
    log_rho_table.DeleteAthenaArray();

  }


  return;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{

  AllocateUserOutputVariables(n_user_var);
  
  // get bottom boundary values, to be used to set up bottom bdr conditions
  Real rbottom = pcoord->x1v(is-1);
  if(rbottom < pmy_mesh->mesh_size.x1min){
    for(int i=1; i<=NGHOST; ++i){
      Real radius = pcoord->x1v(is-i); // note that this means the ghost cell labelled by is-i should have the value bd_data(i-1), see the next comment
      int lleft=0;

      int lright=1;
      while((radius > ini_profile(lright,0)) && (lright < num_of_points_init-1)){
        lright = lright+1;
      }
      if(lright - lleft > 1) lleft = lright -1;

      Real rho = ini_profile(lleft,2) + (radius - ini_profile(lleft,0)) *
                                (ini_profile(lright,2) - ini_profile(lleft,2))
                               /(ini_profile(lright,0) - ini_profile(lleft,0));
      Real tem = ini_profile(lleft,1) + (radius - ini_profile(lleft,0)) *
                                (ini_profile(lright,1) - ini_profile(lleft,1))
                               /(ini_profile(lright,0) - ini_profile(lleft,0));

      Real fr = ini_profile(lleft,5) + (radius - ini_profile(lleft,0)) *
                                (ini_profile(lright,5) - ini_profile(lleft,5))
                               /(ini_profile(lright,0) - ini_profile(lleft,0));

      bd_data(0,i-1) = rho; // you can see here, it is bd_data(i-1) <-> x1v(is-i)
      bd_data(1,i-1) = tem;
      bd_data(2,i-1) = fr;
    }

  }
 
  
  if(NR_RADIATION_ENABLED || IM_RADIATION_ENABLED){
    
      pnrrad->EnrollOpacityFunction(StarOpacity);

  }else{

  }
  

  return;
}


void MeshBlock::UserWorkInLoop(void)
{
  // nothing needs to be done here
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  
  Real gamma = peos->GetGamma();

  NRRadiation *prad = pnrrad;


  //initialize random number
  std::srand(gid);
 
  Real crat, prat;
  if(NR_RADIATION_ENABLED || IM_RADIATION_ENABLED){
    crat = prad->crat;
    prat = prad->prat;
  }else{
    crat = 9166.13;
    prat = 0.0;
  }
  
  Real amp = 0.0;
  const Real seed_amp = pin->DoesParameterExist("problem","seed_amp") ? pin->GetReal("problem","seed_amp") : 0.05;
  
  int kl=ks, ku=ke;
  if(ku > kl){
    ku += NGHOST;
    kl -= NGHOST;
  }
  int jl=js, ju=je;
  if(ju > jl){
    ju += NGHOST;
    jl -= NGHOST;
  }
  int il = is-NGHOST, iu=ie+NGHOST;

  Real rmax = ini_profile(num_of_points_init-1,0);
  Real tem_rmax = ini_profile(num_of_points_init-1,1);
  Real rho_rmax = ini_profile(num_of_points_init-1,2);
  Real mass_rmax = ini_profile(num_of_points_init-1,3);

  Real grav_rmax = grav_const_dimensionless * mass_rmax / pow(rmax,2.0); // all quantities in code units

  // initialize hydro variable
  for(int i=is; i<=ie; ++i) {
    Real &x1 = pcoord->x1v(i); 
    Real tem = tem_rmax;
    Real rho = rho_rmax;   

    // get the position

    Real radflx = 0.0;
    Real ini_profile_max_r = ini_profile(num_of_points_init-1,0);

    if(x1 > ini_profile_max_r){
      tem = tem_rmax;
      Real grav_local = grav_rmax * pow(x1/ini_profile_max_r,2.0);
      rho = rho_rmax * exp(-grav_local*(x1-ini_profile_max_r)/(tem_rmax));
      rho = std::max(rho,rhofloor);

    }else{
      int lleft=0;

      int lright=1;
      while((x1 > ini_profile(lright,0)) && (lright < num_of_points_init-1)){
         lright = lright+1;
      }
      if(lright - lleft > 1) lleft = lright -1;
      
      rho = ini_profile(lleft,2)+ (x1 - ini_profile(lleft,0)) *
                                  (ini_profile(lright,2) - ini_profile(lleft,2))
                                 /(ini_profile(lright,0) - ini_profile(lleft,0));
      tem = ini_profile(lleft,1) + (x1 - ini_profile(lleft,0)) *
                                  (ini_profile(lright,1) - ini_profile(lleft,1))
                                 /(ini_profile(lright,0) - ini_profile(lleft,0));

      radflx = ini_profile(lleft,5) + (x1 - ini_profile(lleft,0)) *
                                  (ini_profile(lright,5) - ini_profile(lleft,5))
                                 /(ini_profile(lright,0) - ini_profile(lleft,0));

    }


    
    if(rho > 0.01 ) amp = seed_amp;
    else amp = 0.0;
    
    
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        phydro->u(IDN,k,j,i) = rho * (1.0 + amp * ((double)rand()/(double)RAND_MAX-0.5));
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS){

          phydro->u(IEN,k,j,i) = tem * rho/(gamma-1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }
        
        if(NR_RADIATION_ENABLED || IM_RADIATION_ENABLED){
          Real er = tem * tem * tem * tem;
          // geometric dialution
          
          for(int ifr=0; ifr<prad->nfreq; ++ifr){
            Real coefa = 0.0, coefb = 0.0;
            for(int n=0; n<prad->nang; ++n){
              // spherical polar coordinate
              Real &miuz = prad->mu(0,k,j,i,n);
              Real &weight = prad->wmu(n);
              if(miuz > 0.0){
                coefa += weight;
                coefb += (miuz * weight);
              }
            }
            
            for(int n=0; n<prad->nang; ++n){
              Real &miuz = prad->mu(0,k,j,i,n);
            
              if(miuz > 0.0){
                prad->ir(k,j,i,ifr*prad->nang+n) = 0.5 *
                                       (er/coefa + radflx/coefb);
              }else{
                prad->ir(k,j,i,ifr*prad->nang+n) = 0.5 *
                                       (er/coefa - radflx/coefb);
              
              }
            
            }
            
          }
        }
      }
    }
  }

  // Opacity will be set during initialization

  
  return;
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
    // custom outputs, here I use the opacity
    // I comment this out because I don't really need it
    /*
    for(int k=ks; k<=ke; k++) {
      for(int j=js; j<=je; j++) {
        for(int i=is; i<=ie; i++) {
          // read dimensionless density (IDN), and dimensionelss gas temparature (IEN) = pressure (IEN) / density (IDN)
          Real rho_dimensionless  = phydro->w(IDN,k,j,i);
          Real T_gas_dimensionless = std::max(phydro->w(IEN,k,j,i)/rho_dimensionless,tfloor);

          // convert to cgs unit, because the interpolation table is based on cgs unit
          Real rho_cgs = rho_dimensionless * rho_unit;
          Real T_gas_cgs = T_gas_dimensionless * T_unit;

          // define Rosseland mean, Planck mean, electron scattering, and absorption opacities
          Real kappa_ross, kappa_planck, kappa_electron, kappa_absorption;

          // interpolate opacities from opacity table
          interpolate_opacity(rho_cgs, T_gas_cgs, kappa_ross, kappa_planck, kappa_electron);

          // output in code unit (1/(rho_unit*l_unit))
          user_out_var(0,k,j,i) = kappa_ross * rho_unit * l_unit;
          user_out_var(1,k,j,i) = kappa_planck * rho_unit * l_unit;
          user_out_var(2,k,j,i) = kappa_electron * rho_unit * l_unit;
        }
      }
    }
    */
  if(NR_RADIATION_ENABLED || IM_RADIATION_ENABLED){
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          user_out_var(0,k,j,i) = pnrrad->pradintegrator->rad_source(0,k,j,i);
          user_out_var(1,k,j,i) = pnrrad->pradintegrator->rad_source(1,k,j,i);
          user_out_var(2,k,j,i) = pnrrad->pradintegrator->rad_source(2,k,j,i);
          user_out_var(3,k,j,i) = pnrrad->pradintegrator->rad_source(3,k,j,i);
        }
      }
    }
  }
}


void StarOpacity(MeshBlock *pmb, AthenaArray<Real> &prim)
{
  // get num of cells
  NRRadiation *prad = pmb->pnrrad;
  int il = pmb->is; int jl = pmb->js; int kl = pmb->ks;
  int iu = pmb->ie; int ju = pmb->je; int ku = pmb->ke;
  il -= NGHOST;
  iu += NGHOST;
  if(ju > jl){
    jl -= NGHOST;
    ju += NGHOST;
  }
  if(ku > kl){
    kl -= NGHOST;
    ku += NGHOST;
  }

  // get opacity for each cell

  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
  for (int i=il; i<=iu; ++i) {
  for (int ifr=0; ifr<prad->nfreq; ++ifr){
    // read dimensionless density (IDN), and dimensionelss gas temparature (IEN) = pressure (IEN) / density (IDN)
    Real rho_dimensionless  = prim(IDN,k,j,i);
    Real T_gas_dimensionless = std::max(prim(IEN,k,j,i)/rho_dimensionless,tfloor);

    // convert to cgs unit, because the interpolation table is based on cgs unit
    Real rho_cgs = rho_dimensionless * rho_unit;
    Real T_gas_cgs = T_gas_dimensionless * T_unit;

    // define Rosseland mean, Planck mean, electron scattering, and absorption opacities
    Real kappa_ross, kappa_planck, kappa_electron, kappa_absorption;

    // interpolate opacities from opacity table
    interpolate_opacity(rho_cgs, T_gas_cgs, kappa_ross, kappa_planck, kappa_electron);

    // absorption is by definintion, Rosseland minus electron scattering 
    // the exact division between absoprtion and scattering is not important if no Compton scattering
    kappa_absorption = kappa_ross - kappa_electron;

    // give rho*kappa to the code, in code unit (1/l_unit)
    prad->sigma_s(k,j,i,ifr) = kappa_electron * rho_cgs * l_unit; //Rosseland mean electron scattering
    prad->sigma_a(k,j,i,ifr) = kappa_absorption * rho_cgs * l_unit; // Rosseland mean absorption
    prad->sigma_p(k,j,i,ifr) = kappa_planck * rho_cgs * l_unit; // Planck mean
    prad->sigma_pe(k,j,i,ifr) = prad->sigma_p(k,j,i,ifr); // also set Planck mean, see Athena++ doc

  }    


 }}}

}


void interpolate_opacity(const Real rho_cgs, const Real T_gas_cgs, Real &kappa_ross, Real &kappa_planck, Real &kappa_electron)
{
  // the interpolation table is in log space, so convert rho and T to log 
  Real log_rho = log10(rho_cgs);
  Real log_T = log10(T_gas_cgs);

  // define log kappa variables
  Real log_kappa_ross, log_kappa_planck, log_kappa_electron;

  // define the indices of the logT/log_rho grid points that is just below/above the desired point to be interpolated
  // the initial values are set to 0, such that if the desired logT/log_rho is below the interpolation range ...
  // ... the following part to find the indices are all skipped, and the indices will be the first indices in the table ...
  // ... corresponding to the smallest logT/log_rho values
  int index_of_tabulated_log_rho_below = 0;
  int index_of_tabulated_log_rho_above = 0;
  int index_of_tabulated_log_T_below = 0;
  int index_of_tabulated_log_T_above = 0;

  // define the Planck opacity scaling factor is temperature is out of range
  Real planck_scaling = 1.0;

  // find below and above indices for log_rho
  while((log_rho > log_rho_table(index_of_tabulated_log_rho_above)) && (index_of_tabulated_log_rho_above < total_num_log_rho_grid - 1))
  {
    index_of_tabulated_log_rho_below = index_of_tabulated_log_rho_above;
    index_of_tabulated_log_rho_above++;
  }
  // in case log_rho is above the last grid in log_rho_table, use the last grid for both indices
  if(index_of_tabulated_log_rho_above == total_num_log_rho_grid - 1 && (log_rho > log_rho_table(index_of_tabulated_log_rho_above)))
    index_of_tabulated_log_rho_below = index_of_tabulated_log_rho_above;


  // find below and above indices for logT
  while((log_T > log_T_table(index_of_tabulated_log_T_above)) && (index_of_tabulated_log_T_above < total_num_log_T_grid - 1))
  {
    index_of_tabulated_log_T_below = index_of_tabulated_log_T_above;
    index_of_tabulated_log_T_above++;
  }
  // in case log_T is above the last grid in log_T_table, use the last grid for both indices
  // set the Planck opacity scaling to be smaller by the ratio T^-3.5
  if(index_of_tabulated_log_T_above == total_num_log_T_grid - 1 && (log_T > log_T_table(index_of_tabulated_log_T_above))){
    index_of_tabulated_log_T_below = index_of_tabulated_log_T_above;
    planck_scaling = pow(10.0, -3.5*(log_T - log_T_table(index_of_tabulated_log_T_above)));
  }

  
  // find the values from the grid
  Real log_kappa_ross_log_rho_below_log_T_below = log_kappa_ross_table(index_of_tabulated_log_T_below,index_of_tabulated_log_rho_below);
  Real log_kappa_ross_log_rho_below_log_T_above = log_kappa_ross_table(index_of_tabulated_log_T_above,index_of_tabulated_log_rho_below);
  Real log_kappa_ross_log_rho_above_log_T_below = log_kappa_ross_table(index_of_tabulated_log_T_below,index_of_tabulated_log_rho_above);
  Real log_kappa_ross_log_rho_above_log_T_above = log_kappa_ross_table(index_of_tabulated_log_T_above,index_of_tabulated_log_rho_above);

  Real log_kappa_planck_log_rho_below_log_T_below = log_kappa_planck_table(index_of_tabulated_log_T_below,index_of_tabulated_log_rho_below);
  Real log_kappa_planck_log_rho_below_log_T_above = log_kappa_planck_table(index_of_tabulated_log_T_above,index_of_tabulated_log_rho_below);
  Real log_kappa_planck_log_rho_above_log_T_below = log_kappa_planck_table(index_of_tabulated_log_T_below,index_of_tabulated_log_rho_above);
  Real log_kappa_planck_log_rho_above_log_T_above = log_kappa_planck_table(index_of_tabulated_log_T_above,index_of_tabulated_log_rho_above);

  Real log_kappa_electron_log_rho_below_log_T_below = log_kappa_electron_table(index_of_tabulated_log_T_below,index_of_tabulated_log_rho_below);
  Real log_kappa_electron_log_rho_below_log_T_above = log_kappa_electron_table(index_of_tabulated_log_T_above,index_of_tabulated_log_rho_below);
  Real log_kappa_electron_log_rho_above_log_T_below = log_kappa_electron_table(index_of_tabulated_log_T_below,index_of_tabulated_log_rho_above);
  Real log_kappa_electron_log_rho_above_log_T_above = log_kappa_electron_table(index_of_tabulated_log_T_above,index_of_tabulated_log_rho_above);

  // find the grid values
  Real tabulated_log_rho_below = log_rho_table(index_of_tabulated_log_rho_below);
  Real tabulated_log_rho_above = log_rho_table(index_of_tabulated_log_rho_above);
  Real tabulated_log_T_below = log_T_table(index_of_tabulated_log_T_below);
  Real tabulated_log_T_above = log_T_table(index_of_tabulated_log_T_above);

  // double interpolation
  if(index_of_tabulated_log_rho_below == index_of_tabulated_log_rho_above) // log_rho out of range
  {
      if(index_of_tabulated_log_T_below == index_of_tabulated_log_T_above)
      // log_rho out of range and log_T out of range, no interpolation needed
      {
        log_kappa_ross = log_kappa_ross_log_rho_below_log_T_below;
        log_kappa_planck = log_kappa_planck_log_rho_below_log_T_below;
        log_kappa_electron = log_kappa_electron_log_rho_below_log_T_below;
      }else{ 
        // log_rho out of range but log_T not out of range, 1d interpolation
        log_kappa_ross = log_kappa_ross_log_rho_below_log_T_below + (log_kappa_ross_log_rho_below_log_T_above - log_kappa_ross_log_rho_below_log_T_below) *
                                (log_T - tabulated_log_T_below)/(tabulated_log_T_above - tabulated_log_T_below);
        log_kappa_planck = log_kappa_planck_log_rho_below_log_T_below + (log_kappa_planck_log_rho_below_log_T_above - log_kappa_planck_log_rho_below_log_T_below) *
                                (log_T - tabulated_log_T_below)/(tabulated_log_T_above - tabulated_log_T_below);
        log_kappa_electron = log_kappa_electron_log_rho_below_log_T_below + (log_kappa_electron_log_rho_below_log_T_above - log_kappa_electron_log_rho_below_log_T_below) *
                                (log_T - tabulated_log_T_below)/(tabulated_log_T_above - tabulated_log_T_below);
      }
    }else{ // log_rho not out of range
      if(index_of_tabulated_log_T_below == index_of_tabulated_log_T_above)
      // log_T out of range but log_rho not out of range, 1d interpolation
      {
        log_kappa_ross = log_kappa_ross_log_rho_below_log_T_below + (log_kappa_ross_log_rho_above_log_T_below - log_kappa_ross_log_rho_below_log_T_below) *
                                (log_rho - tabulated_log_rho_below)/(tabulated_log_rho_above - tabulated_log_rho_below);
        log_kappa_planck = log_kappa_planck_log_rho_below_log_T_below + (log_kappa_planck_log_rho_above_log_T_below - log_kappa_planck_log_rho_below_log_T_below) *
                                (log_rho - tabulated_log_rho_below)/(tabulated_log_rho_above - tabulated_log_rho_below);
        log_kappa_electron = log_kappa_electron_log_rho_below_log_T_below + (log_kappa_electron_log_rho_above_log_T_below - log_kappa_electron_log_rho_below_log_T_below) *
                                (log_rho - tabulated_log_rho_below)/(tabulated_log_rho_above - tabulated_log_rho_below);
      }else{ // neither out of range, 2d interpolation
        log_kappa_ross = log_kappa_ross_log_rho_below_log_T_below * (tabulated_log_T_above - log_T) * (tabulated_log_rho_above - log_rho)/
                               ((tabulated_log_T_above - tabulated_log_T_below) * (tabulated_log_rho_above - tabulated_log_rho_below))
                       + log_kappa_ross_log_rho_below_log_T_above * (log_T - tabulated_log_T_below) * (tabulated_log_rho_above - log_rho)/
                               ((tabulated_log_T_above - tabulated_log_T_below) * (tabulated_log_rho_above - tabulated_log_rho_below))
                       + log_kappa_ross_log_rho_above_log_T_below * (tabulated_log_T_above - log_T) * (log_rho - tabulated_log_rho_below)/
                               ((tabulated_log_T_above - tabulated_log_T_below) * (tabulated_log_rho_above - tabulated_log_rho_below))
                       + log_kappa_ross_log_rho_above_log_T_above * (log_T - tabulated_log_T_below) * (log_rho - tabulated_log_rho_below)/
                               ((tabulated_log_T_above - tabulated_log_T_below) * (tabulated_log_rho_above - tabulated_log_rho_below));
        log_kappa_planck = log_kappa_planck_log_rho_below_log_T_below * (tabulated_log_T_above - log_T) * (tabulated_log_rho_above - log_rho)/
                               ((tabulated_log_T_above - tabulated_log_T_below) * (tabulated_log_rho_above - tabulated_log_rho_below))
                       + log_kappa_planck_log_rho_below_log_T_above * (log_T - tabulated_log_T_below) * (tabulated_log_rho_above - log_rho)/
                               ((tabulated_log_T_above - tabulated_log_T_below) * (tabulated_log_rho_above - tabulated_log_rho_below))
                       + log_kappa_planck_log_rho_above_log_T_below * (tabulated_log_T_above - log_T) * (log_rho - tabulated_log_rho_below)/
                               ((tabulated_log_T_above - tabulated_log_T_below) * (tabulated_log_rho_above - tabulated_log_rho_below))
                       + log_kappa_planck_log_rho_above_log_T_above * (log_T - tabulated_log_T_below) * (log_rho - tabulated_log_rho_below)/
                               ((tabulated_log_T_above - tabulated_log_T_below) * (tabulated_log_rho_above - tabulated_log_rho_below));
        log_kappa_electron = log_kappa_electron_log_rho_below_log_T_below * (tabulated_log_T_above - log_T) * (tabulated_log_rho_above - log_rho)/
                               ((tabulated_log_T_above - tabulated_log_T_below) * (tabulated_log_rho_above - tabulated_log_rho_below))
                       + log_kappa_electron_log_rho_below_log_T_above * (log_T - tabulated_log_T_below) * (tabulated_log_rho_above - log_rho)/
                               ((tabulated_log_T_above - tabulated_log_T_below) * (tabulated_log_rho_above - tabulated_log_rho_below))
                       + log_kappa_electron_log_rho_above_log_T_below * (tabulated_log_T_above - log_T) * (log_rho - tabulated_log_rho_below)/
                               ((tabulated_log_T_above - tabulated_log_T_below) * (tabulated_log_rho_above - tabulated_log_rho_below))
                       + log_kappa_electron_log_rho_above_log_T_above * (log_T - tabulated_log_T_below) * (log_rho - tabulated_log_rho_below)/
                               ((tabulated_log_T_above - tabulated_log_T_below) * (tabulated_log_rho_above - tabulated_log_rho_below));
      }
    }

    kappa_ross = pow(10.0, log_kappa_ross);
    kappa_planck = pow(10.0, log_kappa_planck) * planck_scaling; // multiply Planck opacity by the scaling factor if out of range
    kappa_electron = pow(10.0, log_kappa_electron);

    return;
}


// This function sets boundary condition for primitive variables

// Construct the states seen by the Riemann solver at the inner x1 face.  Reusing the
// Reconstruction methods here is important: for PLM and PPM these generally differ
// from the primitive variables stored at the neighboring cell centers.
void left_and_right_states_at_boundary_hllc(
      MeshBlock *pmb, const AthenaArray<Real> &prim,
      const AthenaArray<Real> &bcc, int k, int j, int is, int ie,
      AthenaArray<Real> &wl, AthenaArray<Real> &wr, X1FaceStates &states)
{
  // Construct the left and right states in exactly the same way as
  // Hydro::CalculateFluxes() in hydro/calculate_fluxes.cpp.
  const int order = pmb->precon->xorder;
  if (order == 1) {
    pmb->precon->DonorCellX1(k, j, is-1, ie+1, prim, bcc, wl, wr);
  } else if (order == 2) {
    pmb->precon->PiecewiseLinearX1(k, j, is-1, ie+1, prim, bcc, wl, wr);
  } else {
    pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, prim, bcc, wl, wr);
  }

  // SoundSpeed() expects a complete, contiguous hydro state array, so make arrays
  // solely for passing the reconstructed left and right states to that interface.
  Real left_hydro_array[NHYDRO], right_hydro_array[NHYDRO];
  for (int n=0; n<NHYDRO; ++n) {
    left_hydro_array[n] = wl(n,is);
    right_hydro_array[n] = wr(n,is);
  }

  states.rho_l = wl(IDN,is);
  states.vel_l = wl(IVX,is);
  states.pressure_l = wl(IPR,is);
  states.sound_speed_l = pmb->peos->SoundSpeed(left_hydro_array);

  states.rho_r = wr(IDN,is);
  states.vel_r = wr(IVX,is);
  states.pressure_r = wr(IPR,is);
  states.sound_speed_r = pmb->peos->SoundSpeed(right_hydro_array);
}

void Inflow_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  // hydro and magnetic bdr conditions
  // set density, velocities, and pressure in ghost cells
  for (int k=ks; k<=ke; ++k) {   // loop over cells in x3 (phi) direction
    for (int j=js; j<=je; ++j) { // loop over cells in x2 (theta) direction
      for (int i=1; i<=ngh; ++i) {  // loop over inner ghost cell in x1 (r) direction
        prim(IDN,k,j,is-i) = bd_data(0,i-1);                 // set density to be the bottom density saved in bd_data, the index is because bd_data is assigned this way
        prim(IVX,k,j,is-i) = 0;                              // set x1 velocity (v_r) to 0
        prim(IVY,k,j,is-i) = 0;                              // set x2 velocity (v_theta) to 0
        prim(IVZ,k,j,is-i) = 0;                              // set x3 velocity (v_phi) to 0
        prim(IEN,k,j,is-i) = bd_data(0,i-1)*bd_data(1,i-1);  // set gas pressure to be the bottom pressure, is density*temperatue in code units        
      }
    }
  }

  // Zero HLLC mass flux is obtained by setting the contact speed to zero.  Including
  // the pressure-dependent wave-speed corrections, the required left face velocity is
  //
  //   u_l = (p_r - p_l - rho_r*c_r*q_r*u_r)/(rho_l*c_l*q_l).
  //
  // The boundary condition controls a cell-centered ghost velocity, whereas HLLC uses
  // a reconstructed face velocity.  Iterate on the innermost ghost-cell value and
  // reconstruct after every correction so that the face value, rather than merely the
  // cell-center value, approaches the expression above.
  AthenaArray<Real> wl, wr;
  wl.NewAthenaArray(NWAVE, pmb->ncells1);
  wr.NewAthenaArray(NWAVE, pmb->ncells1);
  constexpr int max_iterations = 100;
  constexpr Real relative_tolerance = 1.0e-20;
  bool any_boundary_velocity_not_converged = false;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      bool converged = false;
      for (int iteration=0; iteration<max_iterations; ++iteration) {
        X1FaceStates states;
        left_and_right_states_at_boundary_hllc(
            pmb, prim, pmb->pfield->bcc, k, j, is, ie, wl, wr, states);

        const Real rhoa = 0.5*(states.rho_l + states.rho_r);
        const Real ca = 0.5*(states.sound_speed_l + states.sound_speed_r);
        const Real pmid = 0.5*(states.pressure_l + states.pressure_r
                               + (states.vel_l - states.vel_r)*rhoa*ca);
        const Real gamma = pmb->peos->GetGamma();
        const Real ql = (pmid <= states.pressure_l) ? 1.0 :
            std::sqrt(1.0 + (gamma + 1.0)/(2.0*gamma)
                      *(pmid/states.pressure_l - 1.0));
        const Real qr = (pmid <= states.pressure_r) ? 1.0 :
            std::sqrt(1.0 + (gamma + 1.0)/(2.0*gamma)
                      *(pmid/states.pressure_r - 1.0));
        const Real desired_vel_l =
            (states.pressure_r - states.pressure_l
             - states.rho_r*states.sound_speed_r*qr*states.vel_r)
            /(states.rho_l*states.sound_speed_l*ql);
        const Real contact_numerator =
            states.pressure_l - states.pressure_r
            + states.rho_l*states.sound_speed_l*ql*states.vel_l
            + states.rho_r*states.sound_speed_r*qr*states.vel_r;
        const Real pressure_scale =
            std::max(std::max(std::abs(states.pressure_l),
                              std::abs(states.pressure_r)), TINY_NUMBER);
        if (std::abs(contact_numerator) <= relative_tolerance*pressure_scale) {
          converged = true;
          break;
        }

        // Leave prim at a state whose convergence was actually checked rather than
        // applying an unevaluated correction after the final iteration.
        if (iteration < max_iterations-1) {
          prim(IVX,k,j,is-1) += desired_vel_l - states.vel_l;
        }
      }
      if (!converged) any_boundary_velocity_not_converged = true;
    }
  }
  // For VL2, warn only on the half-step boundary call that prepares the second-stage
  // flux.  This produces at most one warning per timestep for this boundary MeshBlock.
  const Real half_step_time = pmb->pmy_mesh->time + 0.5*dt;
  const Real time_tolerance = 1.0e-12
      *std::max(std::max(std::abs(half_step_time), std::abs(dt)), 1.0);
  if (any_boundary_velocity_not_converged
      && std::abs(time - half_step_time) <= time_tolerance) {
    std::printf("### WARNING: inner_x1 HLLC zero-mass-flux velocity did not converge "
                "for at least one boundary cell in MeshBlock gid=%d at time=%.17e\n",
                pmb->gid, time);
  }
  wl.DeleteAthenaArray();
  wr.DeleteAthenaArray();

   // set magnetic field in inlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){
      for(int i=1; i<=ngh; ++i){
        b.x1f(k,j,is-i) = b.x1f(k,j,is);
      }
    }}

    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je+1; ++j){
      for(int i=1; i<=ngh; ++i){
        b.x2f(k,j,is-i) = b.x2f(k,j,is);
      }
    }}

    for(int k=ks; k<=ke+1; ++k){
    for(int j=js; j<=je; ++j){
      for(int i=1; i<=ngh; ++i){
        b.x3f(k,j,is-i) = b.x3f(k,j,is);
      }
    }}
    
  }


  return;
}

void Outflow_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
      // outer hydro bdr conditions in the X1 direction
      // if inflow (V_x1<0), the ghost cells X1 velocities are 0, while density, horizontal velocities, ...
      // ... and pressure are the same as the outer bdr
      // if outflow (V_x1>0), the ghost cells 3D velocities and pressure are the same as the outer bdr, ...
      // ... the density is scaled as r**-2 as extended from outer bdr
{
  for (int k=ks; k<=ke; ++k) {      // loop over cells in x3 (phi) direction
    for (int j=js; j<=je; ++j) {    // loop over cells in x2 (theta) direction
      for (int i=1; i<=ngh; ++i) {  // loop over outer ghost cells in x1 (r) direction
          Real &x1g = pco->x1v(ie+i);  // radius of outer bound of the ghost cell
          Real &x1 = pco->x1v(ie+i-1); // radius of inner bound of the ghost cell, for the inner ghost cell this is the outer bdr
          if(prim(IVX,k,j,ie) < 0.0){  // if inflow
            prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie);  // density is the same as the outer bdr
            prim(IVX,k,j,ie+i) = 0.0;               // X1 velocity is 0
          }else{                       // if outflow
            prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie+i-1) * x1*x1/(x1g*x1g);   // density smaller by a factor of (r_out/r_in)**2 compared to last cell
            prim(IVX,k,j,ie+i) = prim(IVX,k,j,ie);  // X1 velocity is the same as the outer bdr
          }
          prim(IVY,k,j,ie+i) = prim(IVY,k,j,ie);    // X2 velocity is the same as the outer bdr
          prim(IVZ,k,j,ie+i) = prim(IVZ,k,j,ie);    // X3 velocity is the same as the outer bdr
          prim(IEN,k,j,ie+i) = prim(IPR,k,j,ie);    // pressure is the same as the outer bdr
        
      }
    }
  }
  

  return;
}


void Outflow_rad_X1(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
     const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  // outer radiation bdr conditions in the X1 direction
  // simply keep the outflow intensity (vacuum)  
  for (int k=ks; k<=ke; ++k) {                    // loop over cells in x3 (phi) direction
    for (int j=js; j<=je; ++j) {                  // loop over cells in x2 (theta) direction
      for (int i=1; i<=ngh; ++i) {                // loop over inner ghost cells in x1 (r) direction
        for(int ifr=0; ifr<prad->nfreq; ++ifr){   // loop over frequency groups
          for(int n=0; n<prad->nang; ++n){        // loop over angles
            Real miuz = prad->mu(0,k,j,ie+i,ifr*prad->nang+n);    // get intensity direction
            if(miuz > 0.0){                                       // if direction outwards, positive z (r)
              ir(k,j,ie+i,ifr*prad->nang+n) = ir(k,j,ie+i-1,ifr*prad->nang+n);    // keep the intensity
            }else{                                                // if direction inwards, negative z (r)
              ir(k,j,ie+i,ifr*prad->nang+n) = 0.0;                // no incoming intensity
            }  
          }
        }
      }
    }
  }
  

  return;
}


void Inflow_rad_X1(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
     const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
      // inner radiation bdr conditions in the X1 direction
      // inpose a fix density, temperature, luminosity
{
 
  AthenaArray<Real> ir_cm;
  Real *ir_lab;
  ir_cm.NewAthenaArray(prad->n_fre_ang);
 
  for (int k=ks; k<=ke; ++k) {        // loop over cells in x3 (phi) direction
    for (int j=js; j<=je; ++j) {      // loop over cells in x2 (theta) direction
      for (int i=1; i<=ngh; ++i) {    // loop over inner ghost cells

        Real &x1 = pco->x1v(is-i);    
        Real radflx = bd_data(2,i-1);

        Real tem = bd_data(1,i-1);
        Real er = tem * tem * tem * tem;

        
        for(int ifr=0; ifr<prad->nfreq; ++ifr){   // loop over frequency groups
          Real coefa = 0.0, coefb = 0.0;
          for(int n=0; n<prad->nang; ++n){
            // spherical polar coordinate
            Real &miuz = prad->mu(0,k,j,is-i,n); // cos(theta)
            Real &weight = prad->wmu(n); // solid angle
            if(miuz > 0.0){
              coefa += weight;
              coefb += (miuz * weight);
            }
          }

          for(int n=0; n<prad->nang; ++n){ 
            Real &miuz = prad->mu(0,k,j,is-i,n);
            if(miuz > 0.0){
              ir_cm(n) = 0.5 * (er/coefa + radflx/coefb);
              //printf("radflx = %.5e, ", radflx);
              //printf("er = %.5e, ", er);
              //printf("ir_cm+ = %.5e\n", ir_cm(n));
            }
            else{
              ir_cm(n) = 0.5 * (er/coefa - radflx/coefb);
              //printf("radflx = %.5e, ", radflx);
              //printf("er = %.5e, ", er);
              //printf("ir_cm- = %.5e\n", ir_cm(n));
            }
          }
          Real *mux = &(prad->mu(0,k,j,is-i,0));
          Real *muy = &(prad->mu(1,k,j,is-i,0));
          Real *muz = &(prad->mu(2,k,j,is-i,0));

          ir_lab = &(ir(k,j,is-i,0));
    
          prad->pradintegrator->ComToLab(w(IVX,k,j,is-i),w(IVY,k,j,is-i),w(IVZ,k,j,is-i),mux,muy,muz,ir_cm,ir_lab);
         
        }
          
      }
    }
  }
  
  ir_cm.DeleteAthenaArray();
  return;
}
