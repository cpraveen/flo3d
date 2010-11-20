#include <iostream>
#include <string>
#include "grid.h"
#include "parameter.h"

using namespace std;

// Read grid from file
void Grid::read (GridType grid_type, string grid_file)
{
   if(grid_type == gmsh)
      read_gmsh (grid_file);
   else
   {
      cout << "Unknown grid type specified !!!" << endl;
      abort ();
   }

   preproc ();
   info ();
}

void Grid::info ()
{
   cout << "Number of vertices = " << n_vertex;
   cout << "Number of cells    = " << n_cell;
   cout << "Number of faces    = " << n_face;
}
