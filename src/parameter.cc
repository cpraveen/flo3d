#include "parameter.h"

void Parameter::read ()
{
   grid_type = gmsh;
   grid_file = "cube.msh";

   cfl = 0.8;
   max_iter = 1000;
   final_time = 1.0e20;
   min_residue = 1.0e-6;
   mach_inf = 0.6;
}
