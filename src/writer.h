#ifndef __WRITER_H__
#define __WRITER_H__

#include <vector>
#include <string>
#include "grid.h"
#include "material.h"

class Writer
{
   public:
      Writer (const Grid& grid) 
         : 
         grid (&grid),
         has_cell_primitive (false)
         {};
      void attach_vertex_data (std::vector<double>& data, std::string name);
      void attach_cell_data (std::vector<PrimVar>& data);
      void output_vtk (std::string filename);

   private:

      const Grid* grid;

      std::vector< std::vector<double>* > vertex_data;
      std::vector<std::string> vertex_data_name;

      std::vector<PrimVar>* cell_primitive;
      bool has_cell_primitive;

};

#endif
