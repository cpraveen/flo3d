#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include <vector>
#include <string>
#include <map>
#include <fstream>
#include "vec.h"
#include "material.h"
#include "force.h"
#include "ic.h"

// Coefficients for 3-stage RK scheme of Shu-Osher
static const double a_rk[] = {0.0, 3.0/4.0, 1.0/3.0};
static const double b_rk[] = {1.0, 1.0/4.0, 2.0/3.0};

enum GridType {gmsh};

enum BCType { interior, slip, noslip, farfield, inlet, outlet, pressure };

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

      enum ReconstructScheme { first, second, limited };
      ReconstructScheme reconstruct_scheme;

      Material material;

      std::string grid_file;
      GridType    grid_type;

      InitialCondition initial_condition;

      std::map<int,BCType>  bc_type;
      std::map<int,PrimVar> bc_state;

      std::string  write_format;
      unsigned int write_frequency;
      std::vector<std::string> write_variables;
      bool write_restart;
      bool write_vertex_variables;

      std::vector<ForceData> force_data;

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
