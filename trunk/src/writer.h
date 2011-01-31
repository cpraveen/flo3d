#ifndef __WRITER_H__
#define __WRITER_H__

#include <vector>
#include <string>
#include "grid.h"

class Writer
{
   public:
      Writer (const Grid& grid) : grid (&grid) {};
      void attach_data (std::vector<unsigned int>& data, std::string name);
      void output_vtk (std::string filename);

   private:

      const Grid* grid;

      std::vector< std::vector<unsigned int>* > integer_data;
      std::vector<std::string> integer_data_name;

};

#endif
