#include <iostream>
#include <string>
#include <cstdlib>
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
   abort();
}

void Grid::info ()
{
   cout << "Number of vertices = " << n_vertex << endl;
   cout << "Number of cells    = " << n_cell << endl;
   cout << "Number of faces    = " << n_face << endl;
   cout << "Minimum cell volume= " << min_cell_volume << endl;
   cout << "Maximum cell volume= " << max_cell_volume << endl;
}
