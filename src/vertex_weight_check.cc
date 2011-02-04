#include <vector>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "parameter.h"
#include "grid.h"
#include "writer.h"

using namespace std ;

void Grid::vertex_weight_check ()
{     
   unsigned int v;
   vector<double> weight_flag (n_vertex, 1);
   int nega_vertex_weight = 0 ;
   int vertex_boundary = 0 ;

   for (unsigned int i=0; i< n_cell ; i++)
      for (unsigned int j=0; j<4 ; j++)
      {
         v =  cell[i].vertex[j] ;

         if ( cell[i].weight[j] < 0.0 &&  weight_flag[v] == 1)
         {
            weight_flag[v] = 0;
            nega_vertex_weight += 1 ;
         }
      }


   for(unsigned int i=0; i<n_face; ++i)
      for (unsigned int j=0; j<3; j++)
      {
         v =  face[i].vertex[j] ;
         if (weight_flag[v]==0  && face[i].type != -1 )
            {

            vertex_boundary += 1;
            weight_flag[v] = 2 ;

            }          
      }

   cout << " No. of vertex with negative weight: "
        << nega_vertex_weight << endl ;
   cout << " No. of vertex with negative weight on boundary: "
        << vertex_boundary <<endl;
   cout << " No. of vertex with negative weight in interior: "
        << nega_vertex_weight-vertex_boundary << endl;
   

   Writer writer (*this);
   writer.attach_vertex_data (weight_flag, "weight_flag");
   writer.output_vtk ("weight.vtk");
}

