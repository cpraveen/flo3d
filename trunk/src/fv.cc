#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include "parameter.h"
#include "fv.h"
#include "writer.h"

extern bool restart;

using namespace std;

//------------------------------------------------------------------------------
// Set initial condition
//------------------------------------------------------------------------------
void FiniteVolume::initialize ()
{
   primitive.resize (grid.n_cell);
   conserved_old.resize (grid.n_cell);
   primitive_vertex.resize (grid.n_vertex);
   residual.resize (grid.n_cell);
   dt.resize (grid.n_cell);

   if(restart)
   {
      cout << "Reading restart file flo3d.sol ...\n";
      ifstream fi;
      fi.open ("flo3d.sol");
      assert (fi.is_open());
      for(unsigned int i=0; i<grid.n_cell; ++i)
         fi >> primitive[i].density
            >> primitive[i].velocity.x
            >> primitive[i].velocity.y
            >> primitive[i].velocity.z
            >> primitive[i].pressure;
      fi.close ();
   }
   else
   {
      cout << "Setting initial condition to freestream values\n";
      for(unsigned int i=0; i<grid.n_cell; ++i)
         primitive[i] = param.prim_inf;
   }
}

//------------------------------------------------------------------------------
// Interpolate solution from cell center to vertices
//------------------------------------------------------------------------------
void FiniteVolume::interpolate_vertex ()
{
   for(unsigned int i=0; i<grid.n_vertex; ++i)
      primitive_vertex[i].zero ();

   for(unsigned int i=0; i<grid.n_cell; ++i)
      for(unsigned int j=0; j<4; ++j)
         primitive_vertex[grid.cell[i].vertex[j]] += primitive[i] * grid.cell[i].weight[j];
}

//------------------------------------------------------------------------------
// Reconstruct left and right states
// CURRENTLY FIRST ORDER
//------------------------------------------------------------------------------
void FiniteVolume::reconstruct (const unsigned int& f,
                                bool                has_right,
                                vector<PrimVar>&    state) const
{
   // Average on face
   PrimVar face_avg;
   face_avg.zero ();

   for(unsigned int i=0; i<3; ++i)
      face_avg += primitive_vertex[ grid.face[f].vertex[i] ];
   face_avg *= (1.0/3.0);

   // Left state
   unsigned int vl = grid.face[f].lvertex;
   unsigned int cl = grid.face[f].lcell;
   state[0] = primitive[cl] + ( face_avg - primitive_vertex[vl] ) * 0.25;

   // Right state
   if(has_right)
   {
      unsigned int vr = grid.face[f].rvertex;
      unsigned int cr = grid.face[f].rcell;
      state[1] = primitive[cr] + ( face_avg - primitive_vertex[vr] ) * 0.25;
   }
}

//------------------------------------------------------------------------------
// Compute residual for each cell
//------------------------------------------------------------------------------
void FiniteVolume::compute_residual ()
{
   // Interpolate solution from cell to vertex
   interpolate_vertex ();

   for(unsigned int i=0; i<grid.n_cell; ++i)
      residual[i].zero ();

   unsigned int cl, cr;
   vector<PrimVar> state(2);
   Flux flux;

   // Loop over faces and accumulate flux
   for(unsigned int i=0; i<grid.n_face; ++i)
   {
      if(grid.face[i].type == -1)
      {
         cl = grid.face[i].lcell;
         cr = grid.face[i].rcell;
         reconstruct ( i, true, state );
         param.material.num_flux ( state[0], state[1], grid.face[i].normal, flux );
         residual[cl] += flux;
         residual[cr] -= flux;
      }
      else if(param.bc_type[grid.face[i].type] == slip)
      {
         cl = grid.face[i].lcell;
         reconstruct ( i, false, state );
         param.material.slip_flux ( state[0], grid.face[i].normal, flux );
         residual[cl] += flux;
      }
      else if(param.bc_type[grid.face[i].type] == farfield)
      {
         cl = grid.face[i].lcell;
         reconstruct ( i, false, state );
         param.material.num_flux ( state[0], param.prim_inf, grid.face[i].normal, flux );
         residual[cl] += flux;
      }
      else
      {
         cout << "Unknown face type !!!" << endl;
         abort ();
      }
   }
}

//------------------------------------------------------------------------------
// Compute time step
//------------------------------------------------------------------------------
void FiniteVolume::compute_dt ()
{
   for(unsigned int i=0; i<grid.n_cell; ++i)
      dt[i] = 0.0;

   double gamma = param.material.gamma;
   for(unsigned int i=0; i<grid.n_face; ++i)
   {
      double area = grid.face[i].normal.norm();

      unsigned int cl = grid.face[i].lcell;

      double vel_normal_left = primitive[cl].velocity * grid.face[i].normal;
      double c_left = sqrt( gamma * primitive[cl].pressure / primitive[cl].density );

      dt[cl] += fabs(vel_normal_left) + c_left * area;

      if(grid.face[i].type == -1) // Right cell exists only for interior face
      {
         unsigned int cr = grid.face[i].rcell;
         double vel_normal_right = primitive[cr].velocity * grid.face[i].normal;
         double c_right = sqrt( gamma * primitive[cr].pressure / primitive[cr].density );

         dt[cr] += fabs(vel_normal_right) + c_right * area;
      }
   }

   // Compute global time step
   dt_global = 1.0e20;
   for(unsigned int i=0; i<grid.n_cell; ++i)
   {
      dt[i] = param.cfl * grid.cell[i].volume / dt[i];
      dt_global = min( dt_global, dt[i] );
   }

   // For unsteady simulation, use global time step
   if(param.time_mode == "unsteady")
      for(unsigned int i=0; i<grid.n_cell; ++i)
         dt[i] = dt_global;
}

//------------------------------------------------------------------------------
// Matrix-free LUSGS scheme
//------------------------------------------------------------------------------
void FiniteVolume::lusgs (vector<Flux>& update_cell_change,unsigned int & sweep)
{  
   unsigned int j;
   int f;
   ConVar conserved  ;
   vector<PrimVar> primitive_copy(2);
   unsigned int cl, cr;
   vector<PrimVar> state(2);
   Flux flux, flux1, summation_face;
  
   // Forward sweep and backward sweep without interpolation and reconstruction
   // sweep =0 means forward and sweep=1 means backward
  
   if (sweep ==1)
   sweep=grid.n_cell;
   

   for(unsigned int i=sweep; i<grid.n_cell-sweep; ++i)
   {  
      // initially summation over all faces initialized to zero.
      summation_face.mass_flux=0.0;      
      summation_face.energy_flux=0.0;
      summation_face.momentum_flux.x=0.0;      
      summation_face.momentum_flux.y=0.0;		  
      summation_face.momentum_flux.z=0.0;

      dt[i]= (grid.cell[i].volume / dt_global) + 0.5*dt[i];
      j=0;
      while(grid.cell[i].face[j] !=-1 && j<4) 
      {  f= grid.cell[i].face[j];
         if(grid.face[f].type == -1)
	 {
         cr = grid.face[f].rcell;
	 cl = grid.face[f].lcell;
         primitive_copy[0]=primitive[i];
	 

	 //conserved_old[i] = param.material.prim2con (primitive[i]);
	 //reconstruct ( f, true, state );
	 state[0]=primitive[cl];
	 state[1]=primitive[cr];
	 param.material.num_flux ( state[0], state[1], grid.face[f].normal, flux );
	 flux1=flux;
         conserved   = param.material.prim2con (primitive[cl]);
         primitive[cl] = param.material.con2prim(conserved-update_cell_change[cl]*-1);
	 conserved   = param.material.prim2con (primitive[cr]);
	 primitive[cr] = param.material.con2prim(conserved-update_cell_change[cr]*-1);
	 state[0]=primitive[cl];
	 state[1]=primitive[cr];
	 param.material.num_flux ( state[0], state[1], grid.face[f].normal, flux );
	 double gamma = param.material.gamma;
	 double area = grid.face[f].normal.norm();
	 double vel_normal = primitive_copy[0].velocity * grid.face[i].normal;
	 double c  = sqrt( gamma * primitive_copy[0].pressure / primitive_copy[0].density );
         double lamda  = fabs(vel_normal) + c* area;
	
	 summation_face += (update_cell_change[cr]*lamda-(flux-flux1))*(-0.5);
         }
         j=j+1;
	   
     }
      if (sweep ==0) 
      update_cell_change[i] =(residual[i] + summation_face)*(-1/dt[i]);
      
      else if (sweep == grid.n_cell)
      update_cell_change[i]= update_cell_change[i] - summation_face*(0.5/dt[i]);
     }
     for (unsigned int i=0;i<grid.n_cell;i++)
     {
        conserved=conserved_old[i]- (update_cell_change[i])*-1;
        primitive[i]= param.material.con2prim (conserved);
     }
   
}



//------------------------------------------------------------------------------
// Store old conserved variables for multi-stage RK
//------------------------------------------------------------------------------
void FiniteVolume::store_conserved_old ()
{
   for(unsigned int i=0; i<grid.n_cell; ++i)
      conserved_old[i] = param.material.prim2con (primitive[i]);
}

//------------------------------------------------------------------------------
// Update solution by RK scheme
//------------------------------------------------------------------------------
void FiniteVolume::update_solution (const unsigned int r)
{
   double factor;
   ConVar conserved;

   if(param.time_scheme == "rk1" || param.time_scheme == "rk3")
   {
      for(unsigned int i=0; i<grid.n_cell; ++i)
      {
         factor      = dt[i] / grid.cell[i].volume;
         conserved   = param.material.prim2con (primitive[i]);
         conserved   = conserved_old[i] * a_rk[r] +
                       (conserved - residual[i] * factor) * b_rk[r];
         primitive[i]= param.material.con2prim (conserved);
      }
   }
   else if(param.time_scheme == "lusgs")
   {
      // TBD
   }
}

//------------------------------------------------------------------------------
// Compute L2 norm of mass, momentum and energy residuals
//------------------------------------------------------------------------------
void FiniteVolume::compute_residual_norm (const unsigned int iter)
{
   residual_norm.mass_flux     = 0.0;
   residual_norm.momentum_flux = 0.0;
   residual_norm.energy_flux   = 0.0;

   // Sum of squares for each component
   for(unsigned int i=0; i<grid.n_cell; ++i)
   {
      double volume = grid.cell[i].volume;
      residual_norm.mass_flux       += pow(residual[i].mass_flux       / volume, 2);
      residual_norm.momentum_flux.x += pow(residual[i].momentum_flux.x / volume, 2);
      residual_norm.momentum_flux.y += pow(residual[i].momentum_flux.y / volume, 2);
      residual_norm.momentum_flux.z += pow(residual[i].momentum_flux.z / volume, 2);
      residual_norm.energy_flux     += pow(residual[i].energy_flux     / volume, 2);
   }

   // Take square root and normalize by n_cell
   residual_norm.mass_flux       = sqrt (residual_norm.mass_flux)       / grid.n_cell;
   residual_norm.momentum_flux.x = sqrt (residual_norm.momentum_flux.x) / grid.n_cell;
   residual_norm.momentum_flux.y = sqrt (residual_norm.momentum_flux.y) / grid.n_cell;
   residual_norm.momentum_flux.z = sqrt (residual_norm.momentum_flux.z) / grid.n_cell;
   residual_norm.energy_flux     = sqrt (residual_norm.energy_flux)     / grid.n_cell;

   // Total residual of all components
   residual_norm_total = pow(residual_norm.mass_flux, 2) +
                         residual_norm.momentum_flux.square () +
                         pow(residual_norm.energy_flux, 2);
   residual_norm_total = sqrt (residual_norm_total);

   // Copy residual in first iteration for normalization
   if(iter == 1)
   {
      residual_norm_total0 = residual_norm_total;
      cout << "  Initial residual = " << residual_norm_total0 << endl;
      if(residual_norm_total0 == 0.0)
      {
         cout << "  WARNING: Initial residual is zero !!!\n";
         cout << "  WARNING: Setting it to 1.0\n";
         residual_norm_total0 = 1.0;
      }
   }

   residual_norm_total /= residual_norm_total0;

   // File output
   res  << setw(8) << iter << "  " 
        << scientific
        << setprecision (4) 
        << dt_global << "  " 
        << residual_norm_total << "  "
        << residual_norm.mass_flux << "  "
        << residual_norm.momentum_flux.x << "  "
        << residual_norm.momentum_flux.y << "  "
        << residual_norm.momentum_flux.z << "  "
        << residual_norm.energy_flux
        << endl;

   // Screen output
   cout << setw(8) << iter << "  " 
        << scientific
        << setprecision (4) 
        << dt_global << "  " 
        << residual_norm_total << "  "
        << residual_norm.mass_flux << "  "
        << residual_norm.momentum_flux.x << "  "
        << residual_norm.momentum_flux.y << "  "
        << residual_norm.momentum_flux.z << "  "
        << residual_norm.energy_flux
        << endl;
}

//------------------------------------------------------------------------------
// Save solution to file
//------------------------------------------------------------------------------
void FiniteVolume::output (const unsigned int iter)
{
   Writer writer (grid, param.material);
   writer.attach_cell_data (primitive);
   writer.attach_cell_variables (param.write_variables);

   if(param.write_format == "vtk")
      writer.output_vtk ("flo3d.vtk");
}

//------------------------------------------------------------------------------
// Save solution to file for restart
//------------------------------------------------------------------------------
void FiniteVolume::output_restart ()
{
   Writer writer (grid);
   writer.attach_cell_data (primitive);
   writer.output_restart ();
}
//------------------------------------------------------------------------------
// Perform time marching iterations
//------------------------------------------------------------------------------
void FiniteVolume::solve ()
{
   unsigned int iter = 0;
   double time = 0.0;
   residual_norm_total = 1.0e20;
   res.open ("flo3d.res");

   while (residual_norm_total > param.min_residue &&
          iter < param.max_iter && 
          time < param.final_time)
   {
      store_conserved_old ();
      compute_dt ();
      for(unsigned int r=0; r<param.n_rks; ++r)
      {
         compute_residual ();
         update_solution (r);
      }

      ++iter;
      time += dt_global;

      compute_residual_norm (iter);
      if(iter % param.write_frequency == 0) output (iter);
   }

   res.close ();

   // Save final solution
   output (iter);

   if(param.write_restart) output_restart ();

   
}

//------------------------------------------------------------------------------
// This is where the real work starts
//------------------------------------------------------------------------------
void FiniteVolume::run ()
{
   // Read grid from file
   grid.read (param);

   // Set initial condition
   initialize ();

   // Solve the problem
   solve ();
}
