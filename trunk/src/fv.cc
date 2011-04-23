#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include "parameter.h"
#include "fv.h"
#include "writer.h"

extern bool restart;

using namespace std;

//------------------------------------------------------------------------------
// Set initial condition
//------------------------------------------------------------------------------
void FiniteVolume::initialize ()
{
   primitive.resize (grid.n_cell);
   conserved_old.resize (grid.n_cell);
   primitive_vertex.resize (grid.n_vertex);
   residual.resize (grid.n_cell);
   dt.resize (grid.n_cell);

   // For navier stokes, we need gradient of velocity and temperature
   if(param.material.model == "ns")
   {
      dU.resize (grid.n_vertex);
      dV.resize (grid.n_vertex);
      dW.resize (grid.n_vertex);
      dT.resize (grid.n_vertex);
      grid.vertex_volume.resize (grid.n_vertex, 0.0);

      // Compute volume associated with each vertex
      for(unsigned int i=0; i<grid.n_cell; ++i)
         for(unsigned int j=0; j<4; ++j)
            grid.vertex_volume[grid.cell[i].vertex[j]] += (27.0/64.0) * grid.cell[i].volume;
   }

   // If restart option specified, read previous solution from file
   if(restart)
   {
      cout << "Reading restart file flo3d.sol ...\n";
      ifstream fi;
      fi.open ("flo3d.sol");
      assert (fi.is_open());
      for(unsigned int i=0; i<grid.n_cell; ++i)
         fi >> primitive[i].density
            >> primitive[i].velocity.x
            >> primitive[i].velocity.y
            >> primitive[i].velocity.z
            >> primitive[i].pressure;
      fi.close ();
   }
   else
   {
      cout << "Setting initial condition to freestream values\n";
      for(unsigned int i=0; i<grid.n_cell; ++i)
         primitive[i] = param.prim_inf;
   }
}

//------------------------------------------------------------------------------
// Interpolate solution from cell center to vertices
//------------------------------------------------------------------------------
void FiniteVolume::interpolate_vertex ()
{
   for(unsigned int i=0; i<grid.n_vertex; ++i)
      primitive_vertex[i].zero ();

   for(unsigned int i=0; i<grid.n_cell; ++i)
      for(unsigned int j=0; j<4; ++j)
         primitive_vertex[grid.cell[i].vertex[j]] += primitive[i] * grid.cell[i].weight[j];
}

//------------------------------------------------------------------------------
// Compute derivatives of velocity and temperature at grid vertices
//------------------------------------------------------------------------------
void FiniteVolume::compute_vertex_gradients ()
{
   for(unsigned int i=0; i<grid.n_vertex; ++i)
   {
      dU[i] = 0.0;
      dV[i] = 0.0;
      dW[i] = 0.0;
      dT[i] = 0.0;
   }

   for(unsigned int i=0; i<grid.n_face; ++i)
   {
      unsigned int vl = grid.face[i].lvertex;
      unsigned int cl = grid.face[i].lcell;
      double Tl = param.material.temperature (primitive[cl]);

      dU[vl] += grid.face[i].normal * primitive[cl].velocity.x;
      dV[vl] += grid.face[i].normal * primitive[cl].velocity.y;
      dW[vl] += grid.face[i].normal * primitive[cl].velocity.z;
      dT[vl] += grid.face[i].normal * Tl;

      if(grid.face[i].type == -1)
      {  // Interior face
         unsigned int vr = grid.face[i].rvertex;
         unsigned int cr = grid.face[i].rcell;
         double Tr = param.material.temperature (primitive[cr]);

         dU[vr] -= grid.face[i].normal * primitive[cr].velocity.x;
         dV[vr] -= grid.face[i].normal * primitive[cr].velocity.y;
         dW[vr] -= grid.face[i].normal * primitive[cr].velocity.z;
         dT[vr] -= grid.face[i].normal * Tr;
      }
      else
      {  // Boundary face
         unsigned int n0 = grid.face[i].vertex[0];
         unsigned int n1 = grid.face[i].vertex[1];
         unsigned int n2 = grid.face[i].vertex[2];

         // Average state on face
         PrimVar state = primitive_vertex[n0] +
                         primitive_vertex[n1] +
                         primitive_vertex[n2];
         state *= (1.0/3.0);

         // Apply bc to average state
         int face_type = grid.face[i].type;
         BCType bc_type = param.bc_type[grid.face[i].type];

         if(bc_type == noslip)
         {
            state.velocity = 0.0;
         }
         else if(bc_type == slip)
         {
            Vector unit_normal = grid.face[i].normal / grid.face[i].area;
            state.velocity -= unit_normal * (state.velocity * unit_normal);
         }
         else if(bc_type == inlet)
         {
            state = param.bc_state[face_type];
         }
         else if(bc_type == pressure)
         {
            state.pressure = param.bc_state[face_type].pressure;
         }

         double T = param.material.temperature (state);

         // Add contribution to three vertices of face
         dU[n0] += grid.face[i].normal * state.velocity.x;
         dV[n0] += grid.face[i].normal * state.velocity.y;
         dW[n0] += grid.face[i].normal * state.velocity.z;
         dT[n0] += grid.face[i].normal * T;

         dU[n1] += grid.face[i].normal * state.velocity.x;
         dV[n1] += grid.face[i].normal * state.velocity.y;
         dW[n1] += grid.face[i].normal * state.velocity.z;
         dT[n1] += grid.face[i].normal * T;

         dU[n2] += grid.face[i].normal * state.velocity.x;
         dV[n2] += grid.face[i].normal * state.velocity.y;
         dW[n2] += grid.face[i].normal * state.velocity.z;
         dT[n2] += grid.face[i].normal * T;
      }
   }

   // Divide by vertex volume
   for(unsigned int i=0; i<grid.n_vertex; ++i)
   {
      dU[i] *= (9.0/16.0) / grid.vertex_volume[i];
      dV[i] *= (9.0/16.0) / grid.vertex_volume[i];
      dW[i] *= (9.0/16.0) / grid.vertex_volume[i];
      dT[i] *= (9.0/16.0) / grid.vertex_volume[i];
   }

}

//------------------------------------------------------------------------------
// Reconstruct left and right states
//------------------------------------------------------------------------------
void FiniteVolume::reconstruct (const unsigned int& f,
                                bool                has_right,
                                vector<PrimVar>&    state) const
{
   // Average on face
   PrimVar face_avg;
   face_avg.zero ();

   for(unsigned int i=0; i<3; ++i)
      face_avg += primitive_vertex[ grid.face[f].vertex[i] ];
   face_avg *= (1.0/3.0);

   // Left state
   unsigned int vl = grid.face[f].lvertex;
   unsigned int cl = grid.face[f].lcell;
   state[0] = primitive[cl] + ( face_avg - primitive_vertex[vl] ) * 0.25;

   // Right state
   if(has_right)
   {
      unsigned int vr = grid.face[f].rvertex;
      unsigned int cr = grid.face[f].rcell;
      state[1] = primitive[cr] + ( face_avg - primitive_vertex[vr] ) * 0.25;
   }
}

//------------------------------------------------------------------------------
// Compute residual for each cell
//------------------------------------------------------------------------------
void FiniteVolume::compute_residual ()
{
   // Interpolate solution from cell-center to vertex
   interpolate_vertex ();

   // Set residual vector to zero
   for(unsigned int i=0; i<grid.n_cell; ++i)
      residual[i].zero ();

   unsigned int cl, cr;
   vector<PrimVar> state(2);
   Flux flux;

   // Loop over faces and accumulate flux
   for(unsigned int i=0; i<grid.n_face; ++i)
   {
      int face_type = grid.face[i].type;
      BCType bc_type = param.bc_type[grid.face[i].type];
      
      cl = grid.face[i].lcell;

      if(face_type == -1) // interior face
      {
         cr = grid.face[i].rcell;
         reconstruct ( i, true, state );
         param.material.num_flux ( state[0], state[1], grid.face[i].normal, flux );
         residual[cl] += flux;
         residual[cr] -= flux;
      }
      else if(bc_type == slip || bc_type == noslip)
      {
         reconstruct ( i, false, state );
         param.material.slip_flux ( state[0], grid.face[i].normal, flux );
         residual[cl] += flux;
      }
      else if(bc_type == farfield)
      {
         reconstruct ( i, false, state );
         param.material.num_flux (state[0], 
                                  param.bc_state[face_type], 
                                  grid.face[i].normal, 
                                  flux);
         residual[cl] += flux;
      }
      else if(bc_type == inlet)
      {
         param.material.euler_flux(param.bc_state[face_type], 
                                   grid.face[i].normal,
                                   flux);
         residual[cl] += flux;
      }
      else if(bc_type == outlet)
      {
         reconstruct ( i, false, state );
         param.material.euler_flux(state[0],
                                   grid.face[i].normal,
                                   flux);
         residual[cl] += flux;
      }
      else if(bc_type == pressure)
      {
         reconstruct ( i, false, state );
         state[0].pressure = param.bc_state[face_type].pressure;
         param.material.euler_flux(state[0],
                                   grid.face[i].normal,
                                   flux);
         residual[cl] += flux;
      }
      else
      {
         cout << "Unknown face type !!!" << endl;
         abort ();
      }

   }

   // Viscous fluxes
   if(param.material.model == "ns")
   {
      compute_vertex_gradients ();

      for(unsigned int i=0; i<grid.n_face; ++i)
      {
         int face_type = grid.face[i].type;
         BCType bc_type = param.bc_type[face_type];

         unsigned int n0 = grid.face[i].vertex[0];
         unsigned int n1 = grid.face[i].vertex[1];
         unsigned int n2 = grid.face[i].vertex[2];

         // Average derivatives on face
         Vector dUf = (dU[n0] + dU[n1] + dU[n2]) / 3.0;
         Vector dVf = (dV[n0] + dV[n1] + dV[n2]) / 3.0;
         Vector dWf = (dW[n0] + dW[n1] + dW[n2]) / 3.0;
         Vector dTf = (dT[n0] + dT[n1] + dT[n2]) / 3.0;

         // Average state on face
         PrimVar state = primitive_vertex[n0] + 
                         primitive_vertex[n1] + 
                         primitive_vertex[n2];
         state *= (1.0/3.0);

         // Aply bc
         // On solid walls, adiabatic condition
         if(bc_type == noslip)
         {
            state.velocity = 0.0;

            Vector unit_normal = grid.face[i].normal / grid.face[i].area;
            dTf -= unit_normal * (dTf * unit_normal);
         }
         else if(bc_type == slip)
         {
            Vector unit_normal = grid.face[i].normal / grid.face[i].area;
            state.velocity -= unit_normal * (state.velocity * unit_normal);
         }
         else if(bc_type == inlet)
         {
            state = param.bc_state[face_type];
         }
         else if(bc_type == pressure)
         {
            state.pressure = param.bc_state[face_type].pressure;
         }

         // Compute viscous flux
         Flux flux;
         param.material.viscous_flux (state, 
                                      dUf, 
                                      dVf, 
                                      dWf, 
                                      dTf, 
                                      grid.face[i].normal, 
                                      flux);

         // Add flux to residual
         unsigned int cl = grid.face[i].lcell;
         residual[cl] += flux;

         if(face_type == -1)
         {
            unsigned int cr = grid.face[i].rcell;
            residual[cr] -= flux;
         }
      }
   }
}

//------------------------------------------------------------------------------
// Compute time step
//------------------------------------------------------------------------------
void FiniteVolume::compute_dt ()
{
   for(unsigned int i=0; i<grid.n_cell; ++i)
      dt[i] = 0.0;

   double gamma = param.material.gamma;
   for(unsigned int i=0; i<grid.n_face; ++i)
   {
      double area = grid.face[i].area;

      unsigned int cl = grid.face[i].lcell;

      double vel_normal_left = primitive[cl].velocity * grid.face[i].normal;
      double c_left = sqrt( gamma * primitive[cl].pressure / primitive[cl].density );

      dt[cl] += fabs(vel_normal_left) + c_left * area;

      if(grid.face[i].type == -1) // Right cell exists only for interior face
      {
         unsigned int cr = grid.face[i].rcell;
         double vel_normal_right = primitive[cr].velocity * grid.face[i].normal;
         double c_right = sqrt( gamma * primitive[cr].pressure / primitive[cr].density );

         dt[cr] += fabs(vel_normal_right) + c_right * area;
      }
   }

   // Compute global time step
   dt_global = 1.0e20;

   if( param.time_scheme != "lusgs") 
   {
      for(unsigned int i=0; i<grid.n_cell; ++i)
      {
         dt[i] = param.cfl * grid.cell[i].volume / dt[i];
         dt_global = min( dt_global, dt[i] );
      }
   }
   else
      dt_global = 1.0;

   // For unsteady simulation, use global time step
   if(param.time_mode == "unsteady")
   {
      // Adjust time step so that final time is exactly reached
      if(elapsed_time + dt_global > param.final_time)
         dt_global = param.final_time - elapsed_time;
      for(unsigned int i=0; i<grid.n_cell; ++i)
         dt[i] = dt_global;
   }
}

//------------------------------------------------------------------------------
// Matrix-free LUSGS scheme
//------------------------------------------------------------------------------
void FiniteVolume::lusgs ()
{  
   unsigned int f = 0;
   PrimVar prim;
   int neighbour_cell;
   Flux flux_new, flux_old, summation_face; 
   Vector face_normal;
   const double gamma = param.material.gamma;
   const double prandtl = param.material.prandtl;

   // over relaxation factor; higher value improves stability
   // but slows down convergence; needs to be tuned
   const double omega = 1.5;
   
   // Forward sweep
   for (int i=0; i < grid.n_cell; ++i)
   {   
      // initially summation over all faces initialized to zero.
      summation_face.zero ();
      
      // Diagonal scalar value for LUSGS
      // dt contains convective eigenvalue
      dt[i] = dt[i] / param.cfl + 0.5 * omega * dt[i];
      
      // Loop over neighboring cells
      for (unsigned int j=0; j<4; ++j)
      {
         f = grid.cell[i].face[j] ;
         grid.find_cell_neighbour(f, i, neighbour_cell);

         if(param.material.model == "ns")
         {
            double T = param.material.temperature (primitive[i]);
            double mu = param.material.viscosity (T);
            double area = grid.face[f].area;
            double rho = primitive[i].density; // TODO: take average across face

            dt[i] += gamma * mu * area * area / (grid.cell[i].volume * rho * prandtl);
         } 


         if (grid.face[f].type == -1 && neighbour_cell < i)
         {               
            face_normal = grid.face[f].normal;
	         if (grid.face[f].rcell == i)
	            face_normal *= -1.0;
            
            param.material.euler_flux(primitive[neighbour_cell], 
                                      face_normal,
                                      flux_old);
            
            
            PrimVar prim_avg = (primitive[i] + primitive[neighbour_cell])*0.5;
	         double vel_normal  = prim_avg.velocity * face_normal;
	         double c  = sqrt( gamma * prim_avg.pressure / prim_avg.density );
            double area = grid.face[f].area;
	         double lambda  = omega * (fabs(vel_normal) + c * area); 

            // viscous eigenvalue
            if(param.material.model == "ns")
            {
               Vector dr = grid.cell[i].centroid - grid.cell[neighbour_cell].centroid;
               double T = param.material.temperature (prim_avg);
               double mu = param.material.viscosity (T);
               lambda += area * gamma * mu / (dr.norm() * prim_avg.density * prandtl);
            }
	         
            prim = param.material.con2prim(conserved_old[neighbour_cell]
                                           + residual[neighbour_cell]);
            param.material.euler_flux(prim, face_normal, flux_new);
            summation_face += (residual[neighbour_cell]*lambda
                               - (flux_new - flux_old))*(-0.5);
            
         }
         
      } 
            
      residual[i] += summation_face;
      residual[i] *= (-1.0/dt[i]); 
      // Now residual contains increment of conserved variable
   }
   
   // Backward Sweep
   for(int i=grid.n_cell-1; i>=0; --i)
   {  
      // initially summation over all faces initialized to zero.
      summation_face.zero ();

      // Loop over neighboring cells
      for(unsigned int j=0; j<4; ++j)
      {  
         f = grid.cell[i].face[j] ;
         grid.find_cell_neighbour(f , i, neighbour_cell);

         if (grid.face[f].type == -1 && neighbour_cell > i)
         {	                
            face_normal = grid.face[f].normal;
            if ( grid.face[f].rcell == i)
               face_normal *= -1.0;
            	         
            param.material.euler_flux(primitive[neighbour_cell], 
                                      face_normal,
                                      flux_old);
            
            double area = grid.face[f].area;
            
            PrimVar prim_avg = (primitive[i] + primitive[neighbour_cell])*0.5;
            double vel_normal  = prim_avg.velocity * face_normal;
            double c  = sqrt( gamma * prim_avg.pressure / prim_avg.density );
            double lambda  = omega * (fabs(vel_normal) + c * area);

            if(param.material.model == "ns")
            {
               Vector dr = grid.cell[i].centroid - grid.cell[neighbour_cell].centroid;
               double T = param.material.temperature (prim_avg);
               double mu = param.material.viscosity (T);
               lambda += area * gamma * mu / (dr.norm() * prim_avg.density * prandtl);
            }

            
            prim = param.material.con2prim(conserved_old[neighbour_cell] 
                                           + residual[neighbour_cell]);
            param.material.euler_flux(prim, face_normal, flux_new);
            summation_face += (residual[neighbour_cell]*lambda -
                               (flux_new - flux_old))*(-0.5);
         }
      }
      residual[i] -= summation_face * (1.0/dt[i]);
   } 
   
}



//------------------------------------------------------------------------------
// Store old conserved variables for multi-stage RK
//------------------------------------------------------------------------------
void FiniteVolume::store_conserved_old ()
{
   for(unsigned int i=0; i<grid.n_cell; ++i)
      conserved_old[i] = param.material.prim2con (primitive[i]);
}

//------------------------------------------------------------------------------
// Update solution to new time level
//------------------------------------------------------------------------------
void FiniteVolume::update_solution (const unsigned int r)
{
   double factor;
   ConVar conserved;

   if(param.time_scheme == "rk1" || param.time_scheme == "rk3")
   {
      for(unsigned int i=0; i<grid.n_cell; ++i)
      {
         factor      = dt[i] / grid.cell[i].volume;
         conserved   = param.material.prim2con (primitive[i]);
         conserved   = conserved_old[i] * a_rk[r] +
                       (conserved - residual[i] * factor) * b_rk[r];
         primitive[i]= param.material.con2prim (conserved);
      }
   }
   else if(param.time_scheme == "lusgs")
   { 
      // Forward Sweep and backward sweep
      lusgs();
      
      for (unsigned int i=0; i<grid.n_cell; ++i)
      {
         conserved = conserved_old[i] + residual[i];
         primitive[i] = param.material.con2prim (conserved);
      }
   }
}

//------------------------------------------------------------------------------
// Compute L2 norm of mass, momentum and energy residuals
//------------------------------------------------------------------------------
void FiniteVolume::compute_residual_norm (const unsigned int iter)
{
   residual_norm.mass_flux     = 0.0;
   residual_norm.momentum_flux = 0.0;
   residual_norm.energy_flux   = 0.0;

   // Sum of squares for each component
   for(unsigned int i=0; i<grid.n_cell; ++i)
   {
      double volume = grid.cell[i].volume;
      residual_norm.mass_flux       += pow(residual[i].mass_flux       / volume, 2);
      residual_norm.momentum_flux.x += pow(residual[i].momentum_flux.x / volume, 2);
      residual_norm.momentum_flux.y += pow(residual[i].momentum_flux.y / volume, 2);
      residual_norm.momentum_flux.z += pow(residual[i].momentum_flux.z / volume, 2);
      residual_norm.energy_flux     += pow(residual[i].energy_flux     / volume, 2);
   }

   // Take square root and normalize by n_cell
   residual_norm.mass_flux       = sqrt (residual_norm.mass_flux)       / grid.n_cell;
   residual_norm.momentum_flux.x = sqrt (residual_norm.momentum_flux.x) / grid.n_cell;
   residual_norm.momentum_flux.y = sqrt (residual_norm.momentum_flux.y) / grid.n_cell;
   residual_norm.momentum_flux.z = sqrt (residual_norm.momentum_flux.z) / grid.n_cell;
   residual_norm.energy_flux     = sqrt (residual_norm.energy_flux)     / grid.n_cell;

   // Total residual of all components
   residual_norm_total = pow(residual_norm.mass_flux, 2) +
                         residual_norm.momentum_flux.square () +
                         pow(residual_norm.energy_flux, 2);
   residual_norm_total = sqrt (residual_norm_total);

   // Copy residual in first iteration for normalization
   if(iter == 0)
   {
      residual_norm_total0 = residual_norm_total;
      cout << "  Initial residual = " << residual_norm_total0 << endl;
      if(residual_norm_total0 == 0.0)
      {
         cout << "  WARNING: Initial residual is zero !!!\n";
         cout << "  WARNING: Setting it to 1.0\n";
         residual_norm_total0 = 1.0;
      }
   }

   residual_norm_total /= residual_norm_total0;
}

//------------------------------------------------------------------------------
// Log messages to screen and file
//------------------------------------------------------------------------------
void FiniteVolume::log_messages (const unsigned int iter)
{

   if(param.time_mode == "steady")
   {
      // File output
      res  << setw(8) << iter << "  " 
           << scientific
           << setprecision (4) 
           << dt_global << "  " 
           << residual_norm_total << "  "
           << residual_norm.mass_flux << "  "
           << residual_norm.momentum_flux.x << "  "
           << residual_norm.momentum_flux.y << "  "
           << residual_norm.momentum_flux.z << "  "
           << residual_norm.energy_flux
           << endl;

      // Screen output
      cout << setw(8) << iter << "  " 
           << scientific
           << setprecision (4) 
           << dt_global << "  " 
           << residual_norm_total << "  "
           << residual_norm.mass_flux << "  "
           << residual_norm.momentum_flux.x << "  "
           << residual_norm.momentum_flux.y << "  "
           << residual_norm.momentum_flux.z << "  "
           << residual_norm.energy_flux
           << endl;
   }
   else
   {
      // File output
      res  << setw(8) << iter << "  " 
           << scientific
           << setprecision (4) 
           << dt_global << "  " 
           << elapsed_time 
           << endl;

      // Screen output
      cout << setw(8) << iter << "  " 
           << scientific
           << setprecision (4) 
           << dt_global << "  " 
           << elapsed_time 
           << endl;
   }
}

//------------------------------------------------------------------------------
// Save solution to file for visualization
//------------------------------------------------------------------------------
void FiniteVolume::output (const unsigned int iter)
{
   Writer writer (grid, param.material);
   writer.attach_cell_data (primitive);
   writer.attach_cell_variables (param.write_variables);

   if(param.write_format == "vtk")
      writer.output_vtk ("flo3d.vtk");
}

//------------------------------------------------------------------------------
// Save solution to file for restart
//------------------------------------------------------------------------------
void FiniteVolume::output_restart ()
{
   Writer writer (grid);
   writer.attach_cell_data (primitive);
   writer.output_restart ();
}

//------------------------------------------------------------------------------
// Perform time marching iterations
//------------------------------------------------------------------------------
void FiniteVolume::solve ()
{
   unsigned int iter = 0;
   elapsed_time = 0.0;
   residual_norm_total = 1.0e20;
   res.open ("flo3d.res");

   while (residual_norm_total > param.min_residue &&
          iter < param.max_iter && 
          elapsed_time < param.final_time)
   {
      store_conserved_old ();
      compute_dt ();
      for(unsigned int r=0; r<param.n_rks; ++r)
      {
         compute_residual ();

         if(r == param.n_rks-1)
            compute_residual_norm (iter);
         update_solution (r);
      }

      ++iter;
      elapsed_time += dt_global;
      log_messages (iter);

      compute_forces ();
      if(iter % param.write_frequency == 0) output (iter);
   }

   res.close ();

   // Save final solution
   output (iter);

   if(param.write_restart) output_restart ();

   
}

//------------------------------------------------------------------------------
// This is where the real work starts
//------------------------------------------------------------------------------
void FiniteVolume::run ()
{
   // Read grid from file
   grid.read (param);

   create_force_face_list ();

   // Set initial condition
   initialize ();

   // Solve the problem
   solve ();
}
