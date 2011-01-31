#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "writer.h"

using namespace std;

// Add integer data
void Writer::attach_vertex_data (vector<double>& data, string name)
{
   vertex_data.push_back (&data);
   vertex_data_name.push_back (name);
}

void Writer::attach_cell_data (vector<PrimVar>& data)
{
   cell_primitive = &data;
   has_cell_primitive = true;
}

// Write data to vtk file
void Writer::output_vtk (string filename)
{
   ofstream vtk;
   vtk.open (filename.c_str());

   vtk << "# vtk DataFile Version 3.0" << endl;
   vtk << "flo3d" << endl;
   vtk << "ASCII" << endl;
   vtk << "DATASET UNSTRUCTURED_GRID" << endl;
   vtk << "POINTS  " << grid->n_vertex << "  float" << endl;

   for(unsigned int i=0; i<grid->n_vertex; ++i)
      vtk << grid->vertex[i].x << " " 
          << grid->vertex[i].y << " " 
          << grid->vertex[i].z << endl;

   vtk << "CELLS  " << grid->n_cell << " " << 5 * grid->n_cell << endl;
   for(unsigned int i=0; i<grid->n_cell; ++i)
      vtk << 4 << "  " 
          << grid->cell[i].vertex[0] << " "
          << grid->cell[i].vertex[1] << " "
          << grid->cell[i].vertex[2] << " "
          << grid->cell[i].vertex[3] << endl;

   vtk << "CELL_TYPES  " << grid->n_cell << endl;
   for(unsigned int i=0; i<grid->n_cell; ++i)
      vtk << 10 << endl;

   for(unsigned int d=0; d<vertex_data.size(); ++d)
   {
      if(d==0) vtk << "POINT_DATA  " << grid->n_vertex << endl;
      vtk << "SCALARS  " << vertex_data_name[d] << "  float 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_vertex; ++i)
         vtk << (*vertex_data[d])[i] << endl;
   }

   vtk.close ();
}
