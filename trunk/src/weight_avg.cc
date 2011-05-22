#include <vector>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "parameter.h"
#include "grid.h"

using namespace std;

void Grid::weight_average ()
{     
   const double pw = 1.0;

   cout << "Computing averaging weights ...\n";
   
   vector<double> Rx  (n_vertex, 0.0);
   vector<double> Ry  (n_vertex, 0.0);
   vector<double> Rz  (n_vertex, 0.0);
   vector<double> Ixx (n_vertex, 0.0);
   vector<double> Iyy (n_vertex, 0.0);
   vector<double> Izz (n_vertex, 0.0);
   vector<double> Ixy (n_vertex, 0.0);
   vector<double> Iyz (n_vertex, 0.0);
   vector<double> Izx (n_vertex, 0.0);

   vector<int> bdpoint (n_vertex, 0);
   
   // Mark points which are on boundary
   for(unsigned int i=0; i<n_face; ++i)
      if(face[i].type != -1)
         for(unsigned int j=0; j<3; ++j)
            bdpoint[face[i].vertex[j]] = 1;
   
   // Compute matrix entries and rhs
   for (unsigned int i=0; i< n_cell; i++)
      for (unsigned int j=0; j<4; j++)
      {     
         unsigned int v =  cell[i].vertex[j];
         
         Vector dCV = cell[i].centroid - vertex[v];
         dCV       /= pow(dCV.norm(), pw);
         
         Rx[v]  += dCV.x;
         Ry[v]  += dCV.y;
         Rz[v]  += dCV.z;
         Ixx[v] += dCV.x * dCV.x;
         Iyy[v] += dCV.y * dCV.y;
         Izz[v] += dCV.z * dCV.z;
         Ixy[v] += dCV.x * dCV.y;
         Iyz[v] += dCV.y * dCV.z;
         Izx[v] += dCV.z * dCV.x;
      }
   
   // Solve matrix problem
   for (unsigned int i=0; i< n_vertex; i++)
   {
      double Det = Ixx[i]*(Iyy[i]*Izz[i] - Iyz[i]*Iyz[i]) -  
                   Ixy[i]*(Ixy[i]*Izz[i] - Iyz[i]*Izx[i]) +  
                   Izx[i]*(Ixy[i]*Iyz[i] - Iyy[i]*Izx[i]);
      
      // For boundary point, we dont bother with weights since we are 
      // going to do arithmetic averaging
      // WARNING: Smallness of determinant depends on scale
      if (bdpoint[i]==1 && fabs(Det) < 1.0e-15)
      {
         Rx[i] = Ry[i] = Rz[i] = 0.0;
      }  
      else if (bdpoint[i]==0 && fabs(Det) < 1.0e-15)
      {
         cout << "weight_avg: no solution or many solution\n";
         cout << "vertex = " << i << "  Det = " << Det << endl;
         abort();
      }
      else
      {
         double lambda_x = -Rx[i]*(Iyy[i]*Izz[i] - Iyz[i]*Iyz[i]) + 
                            Ry[i]*(Ixy[i]*Izz[i] - Iyz[i]*Izx[i]) -  
                            Rz[i]*(Ixy[i]*Iyz[i] - Iyy[i]*Izx[i]);
         double lambda_y =  Rx[i]*(Ixy[i]*Izz[i] - Izx[i]*Iyz[i]) - 
                            Ry[i]*(Ixx[i]*Izz[i] - Izx[i]*Izx[i]) +  
                            Rz[i]*(Ixx[i]*Iyz[i] - Ixy[i]*Izx[i]);
         double lambda_z = -Rx[i]*(Ixy[i]*Iyz[i] - Izx[i]*Iyy[i]) + 
                            Ry[i]*(Ixx[i]*Iyz[i] - Izx[i]*Ixy[i]) -  
                            Rz[i]*(Ixx[i]*Iyy[i] - Ixy[i]*Ixy[i]);
         
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
            unsigned int v = face[i].vertex[j];
            Rx[v] = Ry[v] = Rz[v] = 0.0;
         }          
   
   
   // Compute normalized weights
   vector<double> sum_weight(n_vertex, 0.0);
   
   for (unsigned int i=0; i<n_cell; i++)
   {  
      for (unsigned int j=0; j<4; j++)
      {     
         unsigned int v = cell[i].vertex[j];
         Vector dCV     = cell[i].centroid - vertex[v];
         dCV           /= pow(dCV.norm(), pw);
         cell[i].weight[j] = 1.0 + Rx[v] * dCV.x +
                                   Ry[v] * dCV.y +
                                   Rz[v] * dCV.z;
         cell[i].weight[j] /= pow(dCV.norm(), pw);
         
         sum_weight[v] += cell[i].weight[j];
      }
   }
   
   // Divide by sum of weights to normalize
   // Also compute min/max weights
   double min_weight =  1.0e20;
   double max_weight = -1.0e20;
   for (unsigned int i=0; i< n_cell; i++)
      for(unsigned int j=0; j<4; ++j)
      {
         cell[i].weight[j] /= sum_weight[cell[i].vertex[j]];
         min_weight = min ( min_weight, cell[i].weight[j] );
         max_weight = max ( max_weight, cell[i].weight[j] );
      }
   
   cout << "  minimum vertex weight                  : "<< min_weight << endl;
   cout << "  maximum vertex weight                  : "<< max_weight << endl;
}
