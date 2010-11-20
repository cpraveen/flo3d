#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include "grid.h"

using namespace std;

void Grid::read_gmsh (const string grid_file)
{
   unsigned int i, count, elem_type, ntags, tag, n_elem;
   string line;

   ifstream file;
   file.open (grid_file.c_str());

   // Skip some lines
   file >> line;
   file >> line;
   file >> line;
   file >> line;

   // Read vertices
   file >> n_vertex;
   assert (n_vertex > 0);
   vertex.resize (n_vertex);

   for(i=0; i<n_vertex; ++i)
      file >> count
           >> vertex[i].x
           >> vertex[i].y
           >> vertex[i].z;

   file >> line;
   file >> line;

   // Read cells
   file >> n_elem;
   assert (n_elem > 0);

   for(i=0; i<n_elem; ++i)
   {
      file >> count
           >> elem_type
           >> ntags;

      assert( ntags==3 );

      if(elem_type == 2) // Triangular face
      {
         face.resize (n_face+1);
         file >> tag >> tag >> tag;
         file >> face[n_face].vertex[0] 
              >> face[n_face].vertex[1] 
              >> face[n_face].vertex[2];
         ++n_face;
      }
      else if(elem_type == 4) // Tetrahedral cell
      {
         cell.resize (n_cell+1);
         file >> tag >> tag >> tag;
         file >> cell[n_cell].vertex[0] 
              >> cell[n_cell].vertex[1] 
              >> cell[n_cell].vertex[2] 
              >> cell[n_cell].vertex[3];
         ++n_cell;
      }
      else
      {
         cout << "Unknown element type !!!" << endl;
         abort ();
      }
   }

   file.close ();

}
