#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include <string>
#include <map>
#include "vec.h"
#include "material.h"

enum GridType {gmsh};

enum BCType { interior, slip, noslip, farfield };

class Parameter
{
   public:
      unsigned int max_iter;
      double cfl;
      double final_time;
      double min_residue;
      double mach_inf;
      Vec    velocity_inf;
      ConVar con_inf;

      std::string grid_file;
      GridType    grid_type;

      std::map<int,BCType> bc;

      void read ();
};

#endif
