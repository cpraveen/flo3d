#include <cassert>
#include "grid.h"

// Compute cell volumes
void Grid::compute_cell_volume ()
{
   unsigned int v0, v1, v2, v3;

   for(unsigned int i=0; i<n_cell; ++i)
   {
      v0 = cell[i].vertex[0];
      v1 = cell[i].vertex[1];
      v2 = cell[i].vertex[2];
      v3 = cell[i].vertex[3];

      cell[i].volume = (vertex[v0] - vertex[v3]) * 
         ( (vertex[v1] - vertex[v3]) ^ (vertex[v2] - vertex[v3]) );
      cell[i].volume /= 6.0;

      assert ( cell[i].volume > 0.0 );
   }
}

// Compute face normals
void Grid::compute_face_normal ()
{
   unsigned int v0, v1, v2;
   Vec r01, r12;

   for(unsigned int i=0; i<n_face; ++i)
   {
      v0 = face[i].vertex[0];
      v1 = face[i].vertex[1];
      v2 = face[i].vertex[2];

      r01 = vertex[v1] - vertex[v0];
      r12 = vertex[v2] - vertex[v1];

      face[i].normal = r01 ^ r12;
   }

}

// Find cells surrounding a cell
void Grid::find_cell_surr_cell ()
{
   unsigned int i;

   // First put all neighbours to -1
   for(i=0; i<n_cell; ++i)
   {
      cell[i].neighbour[0] = -1;
      cell[i].neighbour[1] = -1;
      cell[i].neighbour[2] = -1;
      cell[i].neighbour[3] = -1;
   }


}

// Preprocess the grid
void Grid::preproc ()
{
   compute_cell_volume ();

   find_cell_surr_cell ();

   // make_faces ();
   compute_face_normal ();
}
