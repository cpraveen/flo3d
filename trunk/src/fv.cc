#include <cassert>
#include <cmath>
#include "parameter.h"
#include "fv.h"

// Set initial condition
void FiniteVolume::initialize ()
{
   solution.resize (grid.n_cell);
   solution_old.resize (grid.n_cell);
   solution_vertex.resize (grid.n_vertex);
   residual.resize (grid.n_cell);

   PrimVar prim_var;
   ConVar  con_var;

   prim_var.density    = 1.0;
   prim_var.velocity.x = 1.0;
   prim_var.velocity.y = 0.0;
   prim_var.velocity.z = 0.0;
   prim_var.pressure   = 1.0/(material.gamma * pow(param.mach_inf,2));

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

// Compute residual for each cell
void FiniteVolume::compute_residual ()
{
   interpolate_vertex ();

   for(unsigned int i=0; i<grid.n_cell; ++i)
      residual[i].zero ();
}

// Update solution by RK scheme
void FiniteVolume::update_solution ()
{
   for(unsigned int i=0; i<grid.n_cell; ++i)
      solution[i] = solution_old[i] - residual[i];
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
