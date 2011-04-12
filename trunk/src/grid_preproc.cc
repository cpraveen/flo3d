#include <iostream>
#include <cmath>
#include <cassert>
#include<fstream>
#include<cstdlib>
#include"grid.h"

extern bool debug;

using namespace std;

//------------------------------------------------------------------------------
// Compute cell Centroid
//------------------------------------------------------------------------------
void Grid::compute_cell_centroid ()
{
   for(unsigned int i=0; i<n_cell; ++i)
   { 
      unsigned int v0, v1, v2, v3;
      v0 = cell[i].vertex[0];
      v1 = cell[i].vertex[1];
      v2 = cell[i].vertex[2];
      v3 = cell[i].vertex[3];   
      cell[i].centroid = ( vertex[v0] + vertex[v1] + vertex[v2] + vertex[v3]) / 4.0;
   }   
}

//------------------------------------------------------------------------------
// Compute cell volumes
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
// Compute face normals
//------------------------------------------------------------------------------
void Grid::compute_face_normal_and_area ()
{
   unsigned int v0, v1, v2 ,lvertex ;
   Vector r01, r12;
   double check_normal_l;

   for(unsigned int i=0; i<n_face; ++i)
   {
      v0 = face[i].vertex[0];
      v1 = face[i].vertex[1];
      v2 = face[i].vertex[2];

      r01 = vertex[v1] - vertex[v0];
      r12 = vertex[v2] - vertex[v1];
      
      face[i].normal = (r01 ^ r12) / 2.0;

      // Check orientation of normal
      lvertex        =  face[i].lvertex ;
      check_normal_l = face[i].normal * ( vertex[lvertex]- vertex[v0]);

      if ( check_normal_l > 0.0 )
      	 face[i].normal *=  -1.0 ;

      face[i].area = face[i].normal.norm();
  }
}

//------------------------------------------------------------------------------
// Add new_face to face list
//------------------------------------------------------------------------------
void Grid::add_face (const Face& new_face)
{
   bool found = false;
   unsigned int n = 0;

   // Take any vertex of this face
   unsigned int v = new_face.vertex[0];

   // Loop over all existing faces of vertex v
   while (n<node_face[v].size() && !found)
   {
      unsigned int f = node_face[v][n];

      if(face[f].lcell==-1) // Boundary face not filled yet
      {
         if(face[f] == new_face)
         {
            face[f].lcell   = new_face.lcell;
            face[f].lvertex = new_face.lvertex;

            found = true;
         }
      }
      else if(face[f].rcell == -1) // Boundary or interior face
      {
         if(face[f] == new_face)
         {
            face[f].rcell   = new_face.lcell;
            face[f].rvertex = new_face.lvertex;

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
      face[n_face].lcell      = new_face.lcell; // left cell
      face[n_face].lvertex    = new_face.lvertex; // left vertex
      face[n_face].rcell      = -1; // right cell to be found
      face[n_face].rvertex    = -1; // right vertex to be found

      // Add this face to its three vertices
      for(unsigned int i=0; i<3; ++i)
      {
         v = new_face.vertex[i];
         node_face[v].push_back (n_face);
      }

      ++n_face;
   }

}

//------------------------------------------------------------------------------
// Create interior faces and connectivity data
//------------------------------------------------------------------------------
void Grid::make_faces ()
{
   cout << "Creating faces ..." << endl;
   unsigned int i;

   node_face.resize (n_vertex);

   // Existing boundary faces
   for(i=0; i<n_face; ++i)
   {
      face[i].lcell   = -1;
      face[i].rcell   = -1;
      face[i].lvertex = -1;
      face[i].rvertex = -1;

      // Add this face to the three vertices
      for(unsigned int j=0; j<3; ++j)
      {
         unsigned int v = face[i].vertex[j];
         node_face[v].push_back(i);
      }
   }

   Face new_face;

   for(i=0; i<n_cell; ++i)
   {
      // first face
      new_face.vertex[0] = cell[i].vertex[0];
      new_face.vertex[1] = cell[i].vertex[2];
      new_face.vertex[2] = cell[i].vertex[1];
      new_face.lcell     = i;
      new_face.lvertex   = cell[i].vertex[3];
      add_face (new_face);

      // second face
      new_face.vertex[0] = cell[i].vertex[1];
      new_face.vertex[1] = cell[i].vertex[2];
      new_face.vertex[2] = cell[i].vertex[3];
      new_face.lcell     = i;
      new_face.lvertex   = cell[i].vertex[0];
      add_face (new_face);

      // third face
      new_face.vertex[0] = cell[i].vertex[2];
      new_face.vertex[1] = cell[i].vertex[0];
      new_face.vertex[2] = cell[i].vertex[3];
      new_face.lcell     = i;
      new_face.lvertex   = cell[i].vertex[1];
      add_face (new_face);

      // fourth face
      new_face.vertex[0] = cell[i].vertex[1];
      new_face.vertex[1] = cell[i].vertex[3];
      new_face.vertex[2] = cell[i].vertex[0];
      new_face.lcell     = i;
      new_face.lvertex   = cell[i].vertex[2];
      add_face (new_face);

   }

   cout << "Checking face data ..." << endl;

   // Now check that face data is complete
   for(i=0; i<n_face; ++i)
   {
      // Check left cell exists
      assert(face[i].lcell != -1);

      if(face[i].type == -1) // Interior face, check right cell
         assert(face[i].rcell != -1);
   }

   // Free memory of node_face since we dont need it any more
   for(i=0; i<n_vertex; ++i)
      node_face[i].resize (0);
   node_face.resize (0);
}

//------------------------------------------------------------------------------
// Find cells surrounding a cell
//------------------------------------------------------------------------------
void Grid::find_cell_faces ()
{
	cout << "Finding faces for each cell ..." << endl;
	
	unsigned int i, j;
	int lcell, rcell;
	
	// First put all faces to -1
	for(i=0; i<n_cell; ++i)
	{
		cell[i].face[0] = -1;
		cell[i].face[1] = -1;
		cell[i].face[2] = -1;
		cell[i].face[3] = -1;
	}
   
	for(i=0; i<n_face; ++i)
	{
		lcell = face[i].lcell;
		j = 0;
		while(cell[lcell].face[j] != -1)
			++j;
		cell[lcell].face[j] = i;
				
      rcell = face[i].rcell;
		if(rcell != -1)
		{ 
			j = 0;
			while(cell[rcell].face[j] != -1)
				++j;
			cell[rcell].face[j] = i;
		}
    }
	
}

//------------------------------------------------------------------------------
// Find cell on other side of given face f
//------------------------------------------------------------------------------
void Grid::find_cell_neighbour( const unsigned int& face_no, 
                                const unsigned int& cell_no, 
                                int&                neighbour_cell_no)
{
   if (face[face_no].lcell == cell_no)
      neighbour_cell_no = face[face_no].rcell;
   else if (face[face_no].rcell == cell_no)
      neighbour_cell_no = face[face_no].lcell;
   else
   {
      cout << "find_cell_neighbour: Fatal error !!!\n";
      abort ();
   }
}   



//------------------------------------------------------------------------------
// Renumbering according to cuthill-mckee algorithm
//------------------------------------------------------------------------------
void Grid::renumber_cell()
{  
	int i, j, current_cell, old_cell, k;
	int neighbour =-1;
	vector<Cell> renumbering;
	vector< unsigned int > new_num, old_num;
	// new_num vector says directly the value of renumbering tag for a old cell number
	// old_num vector says what is the value of old cell number for a given renumbering tag
	// renumbering cell vector is a dummy vector to reshuffle all old cell according to new numbering 
	new_num.resize(n_cell,0);
	old_num.resize(n_cell,0);
	renumbering.resize(n_cell);
	
	// Write initial cell numbering to file
	if(debug)
	{
		ofstream out("cell_number_before.dat");
		
		for(i=0; i<n_cell; ++i)
		{ 
			out << i << " " << i << endl;
			for(j=0; j<4; ++j)
			{    
				find_cell_neighbour (cell[i].face[j], i, neighbour);
				if (neighbour != -1)
					out << i << " " << neighbour << endl;
			}
		}
		out.close();
	}
	
	k=1; // k is the renumbering tag according to the algorithm
	for(i=0; i<n_cell; ++i)
	{ 
		j = 0;
		current_cell = old_num[i] ;
		while(cell[current_cell].face[j] != -1 && j <= 3)
		{     
			find_cell_neighbour (cell[current_cell].face[j], current_cell, neighbour);
			old_cell = neighbour;
			if (old_cell != -1)
			{
				if (new_num[old_cell] == 0 && old_cell != 0)
				{
					old_num[k] = old_cell;
					new_num[old_cell] = k;
					++k;
				}
			}
			++j;
		}
	}
	
	// here coping old cells according to new numbers
	for(i=0; i<n_cell; ++i)
	{
		old_cell = old_num[i];
		renumbering[i] = cell[old_cell];
	}
   
	// shifting all renumbered values back to cell
	for(i=0; i<n_cell; ++i)
	{
		cell[i] = renumbering[i];
	}
	
	// changed cell numbers in faces left and right part
	int lcell, rcell;
	for(i=0; i<n_face; ++i)
	{
		lcell = face[i].lcell ;
		rcell = face[i].rcell ;
		face[i].lcell = new_num[lcell];
		if(rcell != -1)
			face[i].rcell = new_num[rcell];
	}         
	
	// Save new cell numbering to file
	if(debug)
	{
		ofstream out("cell_numbered_after.dat");
		for(i=0; i<n_cell; ++i)
		{
			out << i << " " << i << endl;
			for(j=0; j<4; ++j)
			{    
				find_cell_neighbour (cell[i].face[j], i, neighbour);
				if (neighbour != -1)
					out << i << " " << neighbour << endl;
			}
			
		}
		out.close();
	}
}   


//------------------------------------------------------------------------------
// Preprocess the grid
//------------------------------------------------------------------------------
void Grid::preproc ()
{
   compute_cell_centroid ();
   make_faces ();
   find_cell_faces ();
   renumber_cell ();
   compute_cell_volume ();
   compute_face_normal_and_area ();
   weight_average ();
   vertex_weight_check ();
}
