#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include <vector>
#include <string>
#include <map>
#include <fstream>
#include "vec.h"
#include "material.h"

static const double a_rk[] = {0.0, 3.0/4.0, 1.0/3.0};
static const double b_rk[] = {1.0, 1.0/4.0, 2.0/3.0};

enum GridType {gmsh};

enum BCType { interior, slip, noslip, farfield };

class Parameter
{
   public:
      char* file;
      std::ifstream fin;

      std::string time_mode;
      std::string time_scheme;
      unsigned int n_rks;
      unsigned int max_iter;
      double cfl;
      double final_time;
      double min_residue;

      Material material;

      double mach_inf;
      Vector velocity_inf;
      PrimVar prim_inf;

      std::string grid_file;
      GridType    grid_type;

      std::map<int,BCType> bc_type;
      std::map< int,std::vector<double> > bc_state;

      std::string  write_format;
      unsigned int write_frequency;

      void read ();

   private:
      void read_grid ();
      void read_numeric ();
      void read_material ();
      void read_initial_condition ();
      void read_boundary ();
      void read_integrals ();
      void read_output ();
};

#endif
