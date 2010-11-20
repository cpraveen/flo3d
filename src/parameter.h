#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include <string>

enum GridType {gmsh};

class Parameter
{
   public:
      unsigned int max_iter;
      double cfl;
      double final_time;
      double min_residue;
      double mach_inf;

      std::string grid_file;
      GridType grid_type;

      void read ();
};

#endif
