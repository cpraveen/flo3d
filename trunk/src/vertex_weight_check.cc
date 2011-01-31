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
   vector<unsigned int> weight_flag (n_vertex, 1);
   int nega_vertex_weight = 0 ;
   int vertex_interior = 0 ;

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
         if (weight_flag[v]==0  && face[i].type == -1 && weight_flag[v] != 2)
            {

            vertex_interior += 1;
            weight_flag[v] = 2 ;

            }          
      }

   cout << " The number of vertex having negative weight are : "
        << nega_vertex_weight << endl ;
   cout << " The number of vertex having negative weight on interior faces  are : "
        << vertex_interior <<endl;
   cout << " The number of vertex having negative weight on boundary faces  are : "
        << nega_vertex_weight-vertex_interior << endl;
   

   Writer writer (*this);
   writer.attach_data (weight_flag, "weight_flag");
   writer.output_vtk ("weight.vtk");
}

