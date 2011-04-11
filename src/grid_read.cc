#include <iostream>
#include <string>
#include <cstdlib>
#include <cassert>
#include "grid.h"
#include "parameter.h"

using namespace std;

//------------------------------------------------------------------------------
// Read grid from file
//------------------------------------------------------------------------------
void Grid::read (const Parameter& param)
{
   if(param.grid_type == gmsh)
      read_gmsh (param.grid_file);
   else
   {
      cout << "Unknown grid type specified !!!" << endl;
      abort ();
   }

   // At this stage, we have only boundary faces. We save this number.
   n_boundary_face = n_face;

   check_face_type (param.bc_type);
   preproc ();
   info ();
}

//------------------------------------------------------------------------------
// Print some grid information to screen
//------------------------------------------------------------------------------
void Grid::info ()
{
   cout << "Grid information:\n";
   cout << "  Number of vertices   = " << n_vertex << endl;
   cout << "  Number of cells      = " << n_cell << endl;
   cout << "  Number of faces      = " << n_face << endl;
   cout << "  Number of bdry faces = " << n_boundary_face << endl;
   cout << "  Minimum cell volume  = " << min_cell_volume << endl;
   cout << "  Maximum cell volume  = " << max_cell_volume << endl;
}

//------------------------------------------------------------------------------
// Check that all boundary faces have been assigned a bc type
//------------------------------------------------------------------------------
void Grid::check_face_type (const map<int,BCType>& bc_type)
{
   for(unsigned int i=0; i<face.size(); ++i)
   {
      assert (face[i].type != -1);
      if(bc_type.find(face[i].type) == bc_type.end())
      {
         cout << "check_face_type:\n";
         cout << "   No boundary condition specified for\n";
         cout << "   face = " << i << " whose type = " << face[i].type << endl;
         cout << "   There may be more faces with similar problem.\n";
         abort ();
      }
   }
}
