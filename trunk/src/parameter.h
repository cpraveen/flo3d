#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include <vector>
#include <string>
#include <map>
#include "vec.h"
#include "material.h"

static const double a_rk[] = {0.0, 3.0/4.0, 1.0/3.0};
static const double b_rk[] = {1.0, 1.0/4.0, 2.0/3.0};

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
      Vector velocity_inf;
      PrimVar prim_inf;

      std::string grid_file;
      GridType    grid_type;

      std::map<int,BCType> bc_type;
      std::map< int,std::vector<double> > bc_state;

      unsigned int write_frequency;

      void read ();
};

#endif
