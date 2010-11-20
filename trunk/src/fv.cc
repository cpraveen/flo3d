#include <iostream>
#include <cassert>
#include <cmath>
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
   ConVar  con_var;

   prim_var.density  = 1.0;
   prim_var.velocity = param.velocity_inf;
   prim_var.pressure  = 1.0/(GAMMA * pow(param.mach_inf,2));

   con_var = material.prim2con(prim_var);

   for(unsigned int i=0; i<grid.n_cell; ++i)
      solution[i] = con_var;
}

// Interpolate solution from cell center to vertices
void FiniteVolume::interpolate_vertex ()
{
   for(unsigned int i=0; i<grid.n_cell; ++i)
      solution_vertex[i].zero ();

   for(unsigned int i=0; i<grid.n_cell; ++i)
      for(unsigned int j=0; j<4; ++j)
         solution_vertex[grid.cell[i].vertex[j]] += solution[i] * grid.cell[i].weight[j];
}

// Reconstruct left and right states
// CURRENTLY FIRST ORDER
vector<ConVar> FiniteVolume::reconstruct (const unsigned int vl,
                                          const unsigned int cl,
                                          const unsigned int vr,
                                          const unsigned int cr)
{
   vector<ConVar> state(2);

   state[0] = solution[cl];
   state[1] = solution[cr];

   return state;
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
         vl = grid.face[i].nvertex[0];
         vr = grid.face[i].nvertex[1];
         cl = grid.face[i].ncell[0];
         cr = grid.face[i].ncell[1];
         state = reconstruct ( vl, cl, vr, cr );
         flux  = material.num_flux ( state[0], state[1], grid.face[i].normal );
         residual[cl] += flux;
         residual[cr] -= flux;
      }
      else if(param.bc[grid.face[i].type] == slip)
      {
         vl = grid.face[i].nvertex[0];
         cl = grid.face[i].ncell[0];
         flux = material.slip_flux ( solution[cl], grid.face[i].normal );
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
   }
}

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
