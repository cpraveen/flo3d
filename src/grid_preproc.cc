#include <iostream>
#include <cmath>
#include <cassert>
#include "grid.h"

using namespace std;

// Compute cell volumes
void Grid::compute_cell_volume ()
{
   unsigned int v0, v1, v2, v3;

   min_cell_volume =  1.0e20;
   max_cell_volume = -1.0e20;

   for(unsigned int i=0; i<n_cell; ++i)
   {
      v0 = cell[i].vertex[0];
      v1 = cell[i].vertex[1];
      v2 = cell[i].vertex[2];
      v3 = cell[i].vertex[3];

      cell[i].volume = (vertex[v3] - vertex[v0]) * 
         ( (vertex[v1] - vertex[v0]) ^ (vertex[v2] - vertex[v0]) );
      cell[i].volume /= 6.0;

      assert ( cell[i].volume > 0.0 );

      min_cell_volume = min ( min_cell_volume, cell[i].volume );
      max_cell_volume = max ( max_cell_volume, cell[i].volume );
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

// Add new_face to face list
void Grid::add_face (const Face& new_face)
{
   bool found = false;
   unsigned int n = 0;

   // Loop over all existing faces
   while (!found && n < n_face)
   {
      if(face[n].ncell[0]==-1) // Boundary face not filled yet
      {
         if(face[n] == new_face)
         {
            face[n].ncell[0]   = new_face.ncell[0];
            face[n].nvertex[0] = new_face.nvertex[0];

            found = true;
         }
      }
      else if(face[n].ncell[1] == -1) // Boundary or interior face
      {
         if(face[n] == new_face)
         {
            face[n].ncell[1]   = new_face.ncell[0];
            face[n].nvertex[1] = new_face.nvertex[0];

            found = true;
         }
      }

      ++n;
   }

   if(!found) // This is a new face
   {
      face.resize (n_face+1);
      face[n_face].type = -1; // TODO: NEED TO GIVE NEW TYPE FOR INTERIOR FACES
      face[n_face].vertex [0] = new_face.vertex [0];
      face[n_face].vertex [1] = new_face.vertex [1];
      face[n_face].vertex [2] = new_face.vertex [2];
      face[n_face].ncell  [0] = new_face.ncell  [0]; // left cell
      face[n_face].nvertex[0] = new_face.nvertex[0]; // left vertex
      face[n_face].ncell  [1] = -1; // right cell to be found
      face[n_face].nvertex[1] = -1; // right vertex to be found
      ++n_face;
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

// Create interior faces and connectivity data
void Grid::make_faces ()
{
   cout << "Creating faces ..." << endl;
   unsigned int i;

   for(i=0; i<n_face; ++i)
   {
      face[i].ncell[0]   = -1;
      face[i].ncell[1]   = -1;
      face[i].nvertex[0] = -1;
      face[i].nvertex[1] = -1;
   }

   Face new_face;

   for(i=0; i<n_cell; ++i)
   {
      // first face
      new_face.vertex[0] = cell[i].vertex[0];
      new_face.vertex[1] = cell[i].vertex[2];
      new_face.vertex[2] = cell[i].vertex[1];
      new_face.ncell[0]  = i;
      new_face.nvertex[0]= cell[i].vertex[3];
      add_face (new_face);

      // second face
      new_face.vertex[0] = cell[i].vertex[1];
      new_face.vertex[1] = cell[i].vertex[2];
      new_face.vertex[2] = cell[i].vertex[3];
      new_face.ncell[0]  = i;
      new_face.nvertex[0]= cell[i].vertex[0];
      add_face (new_face);

      // third face
      new_face.vertex[0] = cell[i].vertex[2];
      new_face.vertex[1] = cell[i].vertex[0];
      new_face.vertex[2] = cell[i].vertex[3];
      new_face.ncell[0]  = i;
      new_face.nvertex[0]= cell[i].vertex[1];
      add_face (new_face);

      // fourth face
      new_face.vertex[0] = cell[i].vertex[1];
      new_face.vertex[1] = cell[i].vertex[3];
      new_face.vertex[2] = cell[i].vertex[0];
      new_face.ncell[0]  = i;
      new_face.nvertex[0]= cell[i].vertex[2];
      add_face (new_face);

   }

   cout << "Checking face data ..." << endl;

   // Now check that face data is complete
   for(i=0; i<n_face; ++i)
   {
      // Check left cell exists
      assert(face[i].ncell[0] != -1);

      if(face[i].type == -1) // Interior face, check right cell
         assert(face[i].ncell[1] != -1);
   }
}

// Preprocess the grid
void Grid::preproc ()
{
   make_faces ();
   find_cell_surr_cell ();
   compute_cell_volume ();
   compute_face_normal ();
}
