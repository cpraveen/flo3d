#include <cmath>
#include <map>
#include "parameter.h"

using namespace std;

void Parameter::read ()
{
   grid_type = gmsh;
   grid_file = "bump.msh";

   cfl = 0.8;
   max_iter = 10000;
   final_time = 1.0e20;
   min_residue = 1.0e-6;
   mach_inf = 0.6;
   velocity_inf.x = 1.0;
   velocity_inf.y = 0.0;
   velocity_inf.z = 0.0;

   prim_inf.density  = 1.0;
   prim_inf.velocity = velocity_inf;
   prim_inf.pressure  = 1.0/(GAMMA * pow(mach_inf,2));

   bc_type.insert(pair<int,BCType>(100001, farfield));

   bc_type.insert(pair<int,BCType>(100002, slip));

   bc_type.insert(pair<int,BCType>(100003, slip));
   bc_type.insert(pair<int,BCType>(100004, slip));

   bc_type.insert(pair<int,BCType>(100005, slip));
   bc_type.insert(pair<int,BCType>(100006, slip));
   bc_type.insert(pair<int,BCType>(100007, slip));
   bc_type.insert(pair<int,BCType>(100008, farfield));

   write_frequency = 100;
}
