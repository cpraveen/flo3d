#include <vector>
#include "parameter.h"
#include "grid.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std ;

void Grid::weight_average ()
{     
   cout << "Computing averaging weights ...\n";

      unsigned int v ;
      unsigned int v0, v1, v2, v3;

     
      double Xc,Yc,Zc;

      vector<double> Rx  (n_vertex, 0.0);
      vector<double> Ry  (n_vertex, 0.0);
      vector<double> Rz  (n_vertex, 0.0);
      vector<double> Ixx (n_vertex, 0.0);
      vector<double> Iyy (n_vertex, 0.0);
      vector<double> Izz (n_vertex, 0.0);
      vector<double> Ixy (n_vertex, 0.0);
      vector<double> Iyz (n_vertex, 0.0);
      vector<double> Izx (n_vertex, 0.0);

      double dx, dy, dz;

      // Compute matrix entries and rhs
      for (unsigned int i=0; i< n_cell ; i++)
      {  
         v0 = cell[i].vertex[0];
         v1 = cell[i].vertex[1];
         v2 = cell[i].vertex[2];
         v3 = cell[i].vertex[3];
         
         Xc = ( vertex[v0].x + vertex[v1].x + vertex[v2].x + vertex[v3].x ) / 4.0 ;
         Yc = ( vertex[v0].y + vertex[v1].y + vertex[v2].y + vertex[v3].y ) / 4.0 ;
         Zc = ( vertex[v0].z + vertex[v1].z + vertex[v2].z + vertex[v3].z ) / 4.0 ;

         for (unsigned int j=0; j<4 ; j++)
         {     
                v =  cell[i].vertex[j] ;

                dx = Xc - vertex[v].x ;
                dy = Yc - vertex[v].y ;
                dz = Zc - vertex[v].z ;

                Rx[v] += dx;
                Ry[v] += dy;
                Rz[v] += dz;
                Ixx[v] += dx * dx;
                Iyy[v] += dy * dy;
                Izz[v] += dz * dz;
                Ixy[v] += dx * dy;
                Iyz[v] += dy * dz;
                Izx[v] += dz * dx;
          }

       }


      // Solve matrix problem
      double Det, lambda_x, lambda_y, lambda_z;
      for (unsigned int i=0; i< n_vertex ; i++)
      {

       Det = Ixx[i]*(Iyy[i]*Izz[i] - Iyz[i]*Iyz[i]) -  Ixy[i]*(Ixy[i]*Izz[i] - Iyz[i]*Izx[i]) +  Izx[i]*(Ixy[i]*Iyz[i] - Iyy[i]*Izx[i]) ;

       if (Det == 0.0 )
        {
        cout << " System may have no solution or many solution depending on the Dx , Dy , Dz ";
        abort();
        }  
       else
        {
          lambda_x = -Rx[i]*(Iyy[i]*Izz[i] - Iyz[i]*Iyz[i]) + Ry[i]*(Ixy[i]*Izz[i] - Iyz[i]*Izx[i]) -  Rz[i]*(Ixy[i]*Iyz[i] - Iyy[i]*Izx[i]);
          lambda_y =  Rx[i]*(Ixy[i]*Izz[i] - Izx[i]*Iyz[i]) - Ry[i]*(Ixx[i]*Izz[i] - Izx[i]*Izx[i]) +  Rz[i]*(Ixx[i]*Iyz[i] - Ixy[i]*Izx[i]);
          lambda_z = -Rx[i]*(Ixy[i]*Iyz[i] - Izx[i]*Iyy[i]) + Ry[i]*(Ixx[i]*Iyz[i] - Izx[i]*Ixy[i]) -  Rz[i]*(Ixx[i]*Iyy[i] - Ixy[i]*Ixy[i]);

           Rx[i] = lambda_x / Det;
           Ry[i] = lambda_y / Det;
           Rz[i] = lambda_z / Det;
        }

       }
       
      // For boundary vertices, we do arithmetic averaging
       for(unsigned int i=0; i<n_face; ++i)
          if ( face[i].type != -1 )
            for (unsigned int j=0; j<3; j++)
            {
               v =  face[i].vertex[j] ;
               Rx[v] = Ry[v] = Rz[v] = 0.0;
            }          


       // Compute weights
      double min_weight, max_weight ;
      min_weight =  1.0e20;
      max_weight = -1.0e20;
      vector<double> sum_weight(n_vertex, 0.0);

      for (unsigned int i=0; i< n_cell ; i++)
      {  
         v0 = cell[i].vertex[0];
         v1 = cell[i].vertex[1];
         v2 = cell[i].vertex[2];
         v3 = cell[i].vertex[3];
         
         Xc = ( vertex[v0].x + vertex[v1].x + vertex[v2].x + vertex[v3].x ) / 4.0 ;
         Yc = ( vertex[v0].y + vertex[v1].y + vertex[v2].y + vertex[v3].y ) / 4.0 ;
         Zc = ( vertex[v0].z + vertex[v1].z + vertex[v2].z + vertex[v3].z ) / 4.0 ;

         for (unsigned int j=0; j<4 ; j++)
         {     
             v =  cell[i].vertex[j] ;
             cell[i].weight[j] = 1.0 + Rx[v] * (Xc - vertex[v].x) + 
                                       Ry[v] * (Yc - vertex[v].y) + 
                                       Rz[v] * (Zc - vertex[v].z) ;

              min_weight = min ( min_weight, cell[i].weight[j] );
              max_weight = max ( max_weight, cell[i].weight[j] );

              sum_weight[v] += cell[i].weight[j];
           }
        }

      // Divide by sum of weights to normalize
      for (unsigned int i=0; i< n_cell ; i++)
      {  
         for(unsigned int j=0; j<4; ++j)
            cell[i].weight[j] /= sum_weight[cell[i].vertex[j]];
       }


         cout << " The minimum weight calculate on vertex for averageing is "<< min_weight << endl ;
         cout << " The maximum weight calculate on vertex for averageing is "<< max_weight << endl ;

}
