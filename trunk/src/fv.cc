#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include "parameter.h"
#include "fv.h"
#include "writer.h"

using namespace std;

// Set initial condition
void FiniteVolume::initialize ()
{
   primitive.resize (grid.n_cell);
   conserved_old.resize (grid.n_cell);
   primitive_vertex.resize (grid.n_vertex);
   residual.resize (grid.n_cell);
   dt.resize (grid.n_cell);

   for(unsigned int i=0; i<grid.n_cell; ++i)
      primitive[i] = param.prim_inf;
}

// Interpolate solution from cell center to vertices
void FiniteVolume::interpolate_vertex ()
{
   for(unsigned int i=0; i<grid.n_vertex; ++i)
      primitive_vertex[i].zero ();

   for(unsigned int i=0; i<grid.n_cell; ++i)
      for(unsigned int j=0; j<4; ++j)
         primitive_vertex[grid.cell[i].vertex[j]] += primitive[i] * grid.cell[i].weight[j];
}

// Reconstruct left and right states
// CURRENTLY FIRST ORDER
void FiniteVolume::reconstruct (const unsigned int vl,
                                const unsigned int cl,
                                const unsigned int vr,
                                const unsigned int cr,
                                vector<PrimVar>&   state) const
{
   state[0] = primitive[cl];
   state[1] = primitive[cr];
}

// Reconstruct left state for a boundary face
// CURRENTLY FIRST ORDER
void FiniteVolume::reconstruct (const unsigned int vl,
                                const unsigned int cl,
                                PrimVar&           state) const
{
   state = primitive[cl];
}

// Compute residual for each cell
void FiniteVolume::compute_residual ()
{
   // Interpolate solution from cell to vertex
   interpolate_vertex ();

   for(unsigned int i=0; i<grid.n_cell; ++i)
      residual[i].zero ();

   unsigned int vl, vr, cl, cr;
   vector<PrimVar> state(2);
   Flux flux;

   // Loop over faces and accumulate flux
   for(unsigned int i=0; i<grid.n_face; ++i)
   {
      if(param.bc_type[grid.face[i].type] == interior)
      {
         vl = grid.face[i].lvertex;
         vr = grid.face[i].rvertex;
         cl = grid.face[i].lcell;
         cr = grid.face[i].rcell;
         reconstruct ( vl, cl, vr, cr, state );
         material.num_flux ( state[0], state[1], grid.face[i].normal, flux );
         residual[cl] += flux;
         residual[cr] -= flux;
      }
      else if(param.bc_type[grid.face[i].type] == slip)
      {
         vl = grid.face[i].lvertex;
         cl = grid.face[i].lcell;
         reconstruct ( vl, cl, state[0] );
         material.slip_flux ( state[0], grid.face[i].normal, flux );
         residual[cl] += flux;
      }
      else if(param.bc_type[grid.face[i].type] == farfield)
      {
         vl = grid.face[i].lvertex;
         cl = grid.face[i].lcell;
         reconstruct ( vl, cl, state[0] );
         material.num_flux ( state[0], param.prim_inf, grid.face[i].normal, flux );
         residual[cl] += flux;
      }
      else
      {
         cout << "Unknown face type !!!" << endl;
         abort ();
      }
   }
}

// Compute time step
void FiniteVolume::compute_dt ()
{
   for(unsigned int i=0; i<grid.n_cell; ++i)
      dt[i] = 0.0;

   for(unsigned int i=0; i<grid.n_face; ++i)
   {
      double area = grid.face[i].normal.norm();

      unsigned int cl = grid.face[i].lcell;

      double vel_normal_left = primitive[cl].velocity * grid.face[i].normal;
      double c_left = sqrt( GAMMA * primitive[cl].pressure / primitive[cl].density );

      dt[cl] += fabs(vel_normal_left) + c_left * area;

      if(grid.face[i].type == -1) // Right cell exists only for interior face
      {
         unsigned int cr = grid.face[i].rcell;
         double vel_normal_right = primitive[cr].velocity * grid.face[i].normal;
         double c_right = sqrt( GAMMA * primitive[cr].pressure / primitive[cr].density );

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

}

// Store old conserved variables for multi-stage RK
void FiniteVolume::store_conserved_old ()
{
   for(unsigned int i=0; i<grid.n_cell; ++i)
      conserved_old[i] = material.prim2con (primitive[i]);
}

// Update solution by RK scheme
void FiniteVolume::update_solution ()
{
   double factor;
   ConVar conserved;

   for(unsigned int i=0; i<grid.n_cell; ++i)
   {
      factor      = dt[i] / grid.cell[i].volume;
      conserved   = conserved_old[i] - residual[i] * factor;
      primitive[i]= material.con2prim (conserved);
   }
}

// Compute L2 norm of mass, momentum and energy residuals
void FiniteVolume::compute_residual_norm (Flux& norm)
{
   norm.mass_flux     = 0.0;
   norm.momentum_flux = 0.0;
   norm.energy_flux   = 0.0;

   for(unsigned int i=0; i<grid.n_cell; ++i)
   {
      norm.mass_flux       += pow(residual[i].mass_flux,       2);
      norm.momentum_flux.x += pow(residual[i].momentum_flux.x, 2);
      norm.momentum_flux.y += pow(residual[i].momentum_flux.y, 2);
      norm.momentum_flux.z += pow(residual[i].momentum_flux.z, 2);
      norm.energy_flux     += pow(residual[i].energy_flux,     2);
   }
}

// Perform time marching iterations
void FiniteVolume::solve ()
{
   unsigned int iter = 0;
   double time = 0.0;
   Flux   residual_norm;

   while (iter < param.max_iter && time < param.final_time)
   {
      store_conserved_old ();
      compute_dt ();
      compute_residual ();
      update_solution ();
      compute_residual_norm (residual_norm);

      ++iter;
      time += dt_global;
      cout << iter << "  " 
           << dt_global << "  " 
           << residual_norm.mass_flux << "  "
           << residual_norm.momentum_flux.x << "  "
           << residual_norm.momentum_flux.y << "  "
           << residual_norm.momentum_flux.z << "  "
           << residual_norm.energy_flux
           << endl;
   }
}

// Save solution to file
void FiniteVolume::output ()
{
   Writer writer (grid);
   writer.attach_cell_data (primitive);
   writer.output_vtk ("out.vtk");
}

// This is where the real work starts
void FiniteVolume::run ()
{
   param.read ();
   
   // Read grid from file
   grid.read (param);

   // Set initial condition
   initialize ();

   // Solve the problem
   solve ();

   output ();
}
