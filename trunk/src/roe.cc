#include <algorithm>
#include <cmath>
#include "material.h"

using namespace std;

//------------------------------------------------------------------------------
// Numerical flux function
//------------------------------------------------------------------------------
void Material::roe_flux (const PrimVar& left,
                         const PrimVar& right,
                         const Vector& normal,
                         Flux& flux) const
{
   
   double area = normal.norm();
   Vector unit_normal = normal / area;

   // Enthalpy
   double h_left  = gamma*left.pressure/(left.density*(gamma-1.0))
      + 0.5 * left.velocity.square();
   double h_right = gamma*right.pressure/(right.density*(gamma-1.0))
      + 0.5 * right.velocity.square();

   double rho_left_sqrt = sqrt(left.density);
   double rho_right_sqrt = sqrt(right.density);
   double fact_left = rho_left_sqrt / (rho_left_sqrt + rho_right_sqrt);
   double fact_right = 1.0 - fact_left;

   // Roe average state
   double density  = rho_left_sqrt * rho_right_sqrt;
   Vector velocity = left.velocity  * fact_left + 
                     right.velocity * fact_right;
   double h = h_left * fact_left + h_right * fact_right;

   double vel_normal = velocity * unit_normal;
   double c = sqrt( (gamma-1.0) * (h - 0.5*velocity.square()) );

   double dp = right.pressure - left.pressure;
   double vel_left_normal  = left.velocity  * unit_normal;
   double vel_right_normal = right.velocity * unit_normal;
   double dV = vel_right_normal - vel_left_normal;

   if(vel_normal > 0.0)
   {
      double lambda = vel_normal - c;
      double coeff  = 0.5 * (dp - density * c * dV) / (c * c);
      double factor = min(lambda, 0.0) * coeff;

      // Left flux
      flux.mass_flux = left.density * vel_left_normal;
      flux.momentum_flux = unit_normal * left.pressure +
                           left.velocity * left.density * vel_left_normal;
      flux.energy_flux = h_left * flux.mass_flux;

      // Upwind term
      flux.mass_flux     += factor;
      flux.momentum_flux += (velocity - unit_normal * c) * factor;
      flux.energy_flux   += (h - c * vel_normal) * factor;
   }
   else
   {
      double lambda = vel_normal + c;
      double coeff  = 0.5 * (dp + density * c * dV) / (c * c);
      double factor = max(lambda, 0.0) * coeff;

      // Right flux
      flux.mass_flux = right.density * vel_right_normal;
      flux.momentum_flux = unit_normal * right.pressure +
                           right.velocity * right.density * vel_right_normal;
      flux.energy_flux = h_right * flux.mass_flux;

      // Upwind term
      flux.mass_flux     -= factor;
      flux.momentum_flux -= (velocity + unit_normal * c) * factor;
      flux.energy_flux   -= (h + c * vel_normal) * factor;
   }

   flux *= area;
}
