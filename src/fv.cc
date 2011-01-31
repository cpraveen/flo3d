#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include "parameter.h"
#include "fv.h"

using namespace std;

// Set initial condition
void FiniteVolume::initialize ()
{
   solution.resize (grid.n_cell);
   solution_old.resize (grid.n_cell);
   solution_vertex.resize (grid.n_vertex);
   residual.resize (grid.n_cell);
   dt.resize (grid.n_cell);

   PrimVar prim_var;
   prim_var.density  = 1.0;
   prim_var.velocity = param.velocity_inf;
   prim_var.pressure  = 1.0/(GAMMA * pow(param.mach_inf,2));

   param.con_inf = material.prim2con(prim_var);

   for(unsigned int i=0; i<grid.n_cell; ++i)
      solution[i] = param.con_inf;
}

// Interpolate solution from cell center to vertices
void FiniteVolume::interpolate_vertex ()
{
   for(unsigned int i=0; i<grid.n_vertex; ++i)
      solution_vertex[i].zero ();

   for(unsigned int i=0; i<grid.n_cell; ++i)
      for(unsigned int j=0; j<4; ++j)
         solution_vertex[grid.cell[i].vertex[j]] += solution[i] * grid.cell[i].weight[j];
}

// Reconstruct left and right states
// CURRENTLY FIRST ORDER
void FiniteVolume::reconstruct (const unsigned int vl,
                                const unsigned int cl,
                                const unsigned int vr,
                                const unsigned int cr,
                                vector<ConVar>& state) const
{
   state[0] = solution[cl];
   state[1] = solution[cr];
}

// Compute residual for each cell
void FiniteVolume::compute_residual ()
{
   unsigned int vl, vr, cl, cr;
   vector<ConVar> state(2);
   Flux flux;

   // Interpolate solution from cell to vertex
   interpolate_vertex ();

   for(unsigned int i=0; i<grid.n_cell; ++i)
      residual[i].zero ();

   // Loop over faces and accumulate flux
   for(unsigned int i=0; i<grid.n_face; ++i)
   {
      if(param.bc[grid.face[i].type] == interior)
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
      else if(param.bc[grid.face[i].type] == slip)
      {
         vl = grid.face[i].lvertex;
         cl = grid.face[i].lcell;
         flux = material.slip_flux ( solution[cl], grid.face[i].normal );
         residual[cl] += flux;
      }
      else if(param.bc[grid.face[i].type] == farfield)
      {
         vl = grid.face[i].lvertex;
         cl = grid.face[i].lcell;
         material.num_flux ( solution[cl], param.con_inf, grid.face[i].normal, flux );
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
      PrimVar prim_left = material.con2prim(solution[cl]);
      double vel_normal_left = prim_left.velocity * grid.face[i].normal;
      double c_left = sqrt( GAMMA * prim_left.pressure / prim_left.density );

      dt[cl] += fabs(vel_normal_left) + c_left * area;

      if(grid.face[i].type == -1) // Right cell exists only for interior face
      {
         unsigned int cr = grid.face[i].rcell;
         PrimVar prim_right = material.con2prim(solution[cr]);
         double vel_normal_right = prim_right.velocity * grid.face[i].normal;
         double c_right = sqrt( GAMMA * prim_right.pressure / prim_right.density );

         dt[cr] += fabs(vel_normal_right) + c_right * area;
      }
   }

   dt_global = 1.0e20;
   for(unsigned int i=0; i<grid.n_cell; ++i)
   {
      dt[i] = param.cfl * grid.cell[i].volume / dt[i];
      dt_global = min( dt_global, dt[i] );
   }

}

// Update solution by RK scheme
void FiniteVolume::update_solution ()
{
   double factor;

   for(unsigned int i=0; i<grid.n_cell; ++i)
   {
      factor      = dt[i] / grid.cell[i].volume;
      solution[i] = solution_old[i] - residual[i] * factor;
   }
}

// Perform time marching iterations
void FiniteVolume::solve ()
{
   unsigned int iter = 0;
   double time = 0.0;

   while (iter < param.max_iter && time < param.final_time)
   {
      solution_old = solution;

      compute_residual ();

      update_solution ();

      ++iter;
      time += dt_global;
   }
}

// Save solution to file
void FiniteVolume::output ()
{
}

// This is where the real work starts
void FiniteVolume::run ()
{
   param.read ();
   
   // Read grid from file
   grid.read (param.grid_type, param.grid_file);

   // Set initial condition
   initialize ();

   // Solve the problem
   solve ();

   output ();
}
