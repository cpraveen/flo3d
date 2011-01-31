#include <vector>
#include "parameter.h"
#include "grid.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std ;

void Grid::vertex_weight_check ()
{     
   unsigned int i,j, Val;
   vector<int> weight_flag (n_vertex, 1);
   int nega_vertex_weight = 0 ;
   int vertex_interior = 0 ;

   for ( i=0; i< n_cell ; i++)
      for ( j=0; j<4 ; j++)
      {
         Val =  cell[i].vertex[j] ;

         if ( cell[i].weight[j] < 0.0 &&  weight_flag[Val] == 1)
         {
            weight_flag[Val] = 0;
            nega_vertex_weight += 1 ;
         }
      }


   for(i=0; i<n_face; ++i)
      for ( j=0; j<3 ; j++)
      {
         Val =  face[i].vertex[j] ;
         if (weight_flag[Val]==0  && face[i].type == -1 && weight_flag[Val] != 2)
            {

            vertex_interior += 1;
            weight_flag[Val] = 2 ;

            }          
      }

   cout << " The number of vertex having negative weight are : "
        << nega_vertex_weight << endl ;
   cout << " The number of vertex having negative weight on interior faces  are : "
        << vertex_interior <<endl;
   cout << " The number of vertex having negative weight on boundary faces  are : "
        << nega_vertex_weight-vertex_interior << endl;
   
}

