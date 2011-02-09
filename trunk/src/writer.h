#ifndef __WRITER_H__
#define __WRITER_H__

#include <vector>
#include <string>
#include "grid.h"
#include "material.h"

class Writer
{
   public:
      Writer (const Grid&     grid)
         : 
         grid (&grid),
         has_vertex_primitive (false),
         has_cell_primitive (false),
         write_cell_density (false),
         write_cell_velocity (false),
         write_cell_pressure (false),
         write_cell_mach (false)
         {};
      Writer (const Grid&     grid,
              const Material& material) 
         : 
         grid (&grid),
         material (&material),
         has_vertex_primitive (false),
         has_cell_primitive (false),
         write_cell_density (false),
         write_cell_velocity (false),
         write_cell_pressure (false),
         write_cell_mach (false)
         {};
      void attach_vertex_data (std::vector<double>& data, std::string name);
      void attach_vertex_data (std::vector<PrimVar>& data);
      void attach_cell_data (std::vector<PrimVar>& data);
      void attach_cell_variables (const std::vector<std::string>& variables);
      void output_vtk (std::string filename);
      void output_restart ();

   private:

      const Grid*     grid;
      const Material* material;

      std::vector< std::vector<double>* > vertex_data;
      std::vector<std::string> vertex_data_name;
      std::vector<PrimVar>* vertex_primitive;
      bool has_vertex_primitive;

      std::vector<PrimVar>* cell_primitive;
      bool has_cell_primitive;

      bool write_cell_density;
      bool write_cell_velocity;
      bool write_cell_pressure;
      bool write_cell_mach;

};

#endif
