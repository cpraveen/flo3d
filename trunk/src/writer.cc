#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include "writer.h"

using namespace std;

// Add integer data
void Writer::attach_vertex_data (vector<double>& data, string name)
{
   vertex_data.push_back (&data);
   vertex_data_name.push_back (name);
}

// Add primitive variables defined at vertices
void Writer::attach_vertex_data (vector<PrimVar>& data)
{
   assert (!has_vertex_primitive);
   vertex_primitive = &data;
   has_vertex_primitive = true;
}

// Add primitive variables defined at cells
void Writer::attach_cell_data (vector<PrimVar>& data)
{
   assert (!has_cell_primitive);
   cell_primitive = &data;
   has_cell_primitive = true;
}

// Add mach number at cell centers
void Writer::attach_cell_mach ()
{
   assert (!has_cell_mach);
   assert ( has_cell_primitive);
   has_cell_mach = true;
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

   // Write vertex data
   if(vertex_data.size() > 0 || has_vertex_primitive) 
      vtk << "POINT_DATA  " << grid->n_vertex << endl;

   // Write vertex data to file
   for(unsigned int d=0; d<vertex_data.size(); ++d)
   {
      vtk << "SCALARS  " << vertex_data_name[d] << "  float 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_vertex; ++i)
         vtk << (*vertex_data[d])[i] << endl;
   }

   // If vertex primitive data is available, write to file
   if (has_vertex_primitive)
   {
      vtk << "SCALARS density float 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_cell; ++i)
         vtk << (*vertex_primitive)[i].density << endl;

      vtk << "SCALARS pressure float 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_cell; ++i)
         vtk << (*vertex_primitive)[i].pressure << endl;

      vtk << "VECTORS velocity float" << endl;
      for(unsigned int i=0; i<grid->n_cell; ++i)
         vtk << (*vertex_primitive)[i].velocity.x << "  "
             << (*vertex_primitive)[i].velocity.y << "  "
             << (*vertex_primitive)[i].velocity.z
             << endl;
   }

   // Write cell data
   if(has_cell_primitive || has_cell_mach)
      vtk << "CELL_DATA  " << grid->n_cell << endl;

   // If cell primitive data is available, write to file
   if (has_cell_primitive)
   {
      vtk << "SCALARS density float 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_cell; ++i)
         vtk << (*cell_primitive)[i].density << endl;

      vtk << "SCALARS pressure float 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_cell; ++i)
         vtk << (*cell_primitive)[i].pressure << endl;

      vtk << "VECTORS velocity float" << endl;
      for(unsigned int i=0; i<grid->n_cell; ++i)
         vtk << (*cell_primitive)[i].velocity.x << "  "
             << (*cell_primitive)[i].velocity.y << "  "
             << (*cell_primitive)[i].velocity.z
             << endl;
   }

   if(has_cell_mach)
   {
      vtk << "SCALARS mach float 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_cell; ++i)
      {
         double sonic_square = GAMMA * (*cell_primitive)[i].pressure /
                                       (*cell_primitive)[i].density;
         double mach = sqrt ( (*cell_primitive)[i].velocity. square() /
                               sonic_square );
         vtk << mach << endl;
      }
   }

   vtk.close ();
}
