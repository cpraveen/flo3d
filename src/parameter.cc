#include <cmath>
#include "parameter.h"

using namespace std;

void Parameter::read ()
{
   grid_type = gmsh;
   grid_file = "cube.msh";

   cfl = 0.8;
   max_iter = 1000;
   final_time = 1.0e20;
   min_residue = 1.0e-6;
   mach_inf = 0.6;
   velocity_inf.x = 1.0;
   velocity_inf.y = 0.0;
   velocity_inf.z = 0.0;

   prim_inf.density  = 1.0;
   prim_inf.velocity = velocity_inf;
   prim_inf.pressure  = 1.0/(GAMMA * pow(mach_inf,2));

   bc_type.insert(pair<int,BCType>(10001, slip));
}
