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

      unsigned int i,j, Val ;
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

      for ( i=0; i< n_cell ; i++)
      {  
         v0 = cell[i].vertex[0];
         v1 = cell[i].vertex[1];
         v2 = cell[i].vertex[2];
         v3 = cell[i].vertex[3];
         
         Xc = ( vertex[v0].x + vertex[v1].x + vertex[v2].x + vertex[v3].x ) / 4.0 ;
         Yc = ( vertex[v0].y + vertex[v1].y + vertex[v2].y + vertex[v3].y ) / 4.0 ;
         Zc = ( vertex[v0].z + vertex[v1].z + vertex[v2].z + vertex[v3].z ) / 4.0 ;

         for ( j=0; j<4 ; j++)
         {     
                Val =  cell[i].vertex[j] ;

                dx = Xc - vertex[Val].x ;
                dy = Yc - vertex[Val].y ;
                dz = Zc - vertex[Val].z ;

                Rx[Val] += dx;
                Ry[Val] += dy;
                Rz[Val] += dz;
                Ixx[Val] += dx * dx;
                Iyy[Val] += dy * dy;
                Izz[Val] += dz * dz;
                Ixy[Val] += dx * dy;
                Iyz[Val] += dy * dz;
                Izx[Val] += dz * dx;
          }

       }


      double Det, lambda_x, lambda_y, lambda_z;
      for ( i=0; i< n_vertex ; i++)
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
       
       unsigned int v;

       for(unsigned int i=0; i<n_face; ++i)
          if ( face[i].type != -1 )
            for (unsigned int j=0; j<3; j++)
            {
               v =  face[i].vertex[j] ;
               Rx[v] = 0.0;
               Ry[v] = 0.0;
               Rz[v] = 0.0;
            }          


      double min_weight , max_weight ;
      min_weight =  1.0e20;
      max_weight = -1.0e20;
      vector<double> sum_weight(n_vertex, 0.0);

      for ( i=0; i< n_cell ; i++)
      {  
         v0 = cell[i].vertex[0];
         v1 = cell[i].vertex[1];
         v2 = cell[i].vertex[2];
         v3 = cell[i].vertex[3];
         
         Xc = ( vertex[v0].x + vertex[v1].x + vertex[v2].x + vertex[v3].x ) / 4.0 ;
         Yc = ( vertex[v0].y + vertex[v1].y + vertex[v2].y + vertex[v3].y ) / 4.0 ;
         Zc = ( vertex[v0].z + vertex[v1].z + vertex[v2].z + vertex[v3].z ) / 4.0 ;

         for ( j=0; j<4 ; j++)
         {     
             Val =  cell[i].vertex[j] ;
             cell[i].weight[j] = 1.0 + Rx[Val] * (Xc - vertex[Val].x) + 
                                       Ry[Val] * (Yc - vertex[Val].y) + 
                                       Rz[Val] * (Zc - vertex[Val].z) ;

              min_weight = min ( min_weight, cell[i].weight[j] );
              max_weight = max ( max_weight, cell[i].weight[j] );

              sum_weight[Val] += cell[i].weight[j];
           }
        }

      // Divide by sum of weights to normalize
      for ( i=0; i< n_cell ; i++)
      {  
         for(j=0; j<4; ++j)
            cell[i].weight[j] /= sum_weight[cell[i].vertex[j]];
       }


         cout << " The minimum weight calculate on vertex for averageing is "<< min_weight << endl ;
         cout << " The maximum weight calculate on vertex for averageing is "<< max_weight << endl ;

}
