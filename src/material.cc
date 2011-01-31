#include <algorithm>
#include <cmath>
#include "material.h"

using namespace std;

ConVar& ConVar::operator= (const ConVar& con_var)
{
   density  = con_var.density;
   momentum = con_var.momentum;
   energy   = con_var.energy;

   return *this;
}

ConVar& ConVar::operator+= (const ConVar& con_var)
{
   density  += con_var.density;
   momentum += con_var.momentum;
   energy   += con_var.energy;

   return *this;
}

ConVar ConVar::operator+ (const ConVar& con_var) const
{
   ConVar result;

   result.density  = density  + con_var.density;
   result.momentum = momentum + con_var.momentum;
   result.energy   = energy   + con_var.energy;

   return result;
}

ConVar ConVar::operator- (const ConVar& con_var) const
{
   ConVar result;

   result.density  = density  - con_var.density;
   result.momentum = momentum - con_var.momentum;
   result.energy   = energy   - con_var.energy;

   return result;
}

ConVar ConVar::operator- (const Flux& flux) const
{
   ConVar result;

   result.density  = density  - flux.mass_flux;
   result.momentum = momentum - flux.momentum_flux;
   result.energy   = energy   - flux.energy_flux;

   return result;
}

ConVar ConVar::operator* (const double scalar) const
{
   ConVar result;

   result.density  = density  * scalar;
   result.momentum = momentum * scalar; 
   result.energy   = energy   * scalar;

   return result;
}

PrimVar PrimVar::operator* (const double scalar) const
{
   PrimVar result;

   result.density  = density  * scalar;
   result.velocity = velocity * scalar; 
   result.pressure = pressure * scalar;

   return result;
}

PrimVar& PrimVar::operator+= (const PrimVar& prim_var)
{
   density  += prim_var.density;
   velocity += prim_var.velocity;
   pressure += prim_var.pressure;

   return *this;
}

// Set conserved variable to zero
void PrimVar::zero ()
{
   density  = 0.0;
   velocity = 0.0;
   pressure = 0.0;
}

double ConVar::pressure () const
{
   return (GAMMA - 1.0) * (energy - 0.5 * momentum.square() / density);
}

void Flux::zero ()
{
   mass_flux     = 0.0;
   momentum_flux = 0.0;
   energy_flux   = 0.0;
}

Flux& Flux::operator+= (const Flux& flux)
{
   mass_flux     += flux.mass_flux;
   momentum_flux += flux.momentum_flux;
   energy_flux   += flux.energy_flux;

   return *this;
}

Flux& Flux::operator-= (const Flux& flux)
{
   mass_flux     -= flux.mass_flux;
   momentum_flux -= flux.momentum_flux;
   energy_flux   -= flux.energy_flux;

   return *this;
}

// Multiply flux by a scalar
Flux& Flux::operator*= (const double& scalar)
{
   mass_flux     *= scalar;
   momentum_flux *= scalar;
   energy_flux   *= scalar;

   return *this;
}

// Add two fluxes
Flux Flux::operator+ (const Flux& flux)
{
   Flux result;

   result.mass_flux     = mass_flux + flux.mass_flux;
   result.momentum_flux = momentum_flux + flux.momentum_flux;
   result.energy_flux   = energy_flux + flux.energy_flux;

   return result;
}

Flux Flux::operator* (const double scalar)
{
   Flux result;

   result.mass_flux     = mass_flux * scalar;
   result.momentum_flux = momentum_flux * scalar;
   result.energy_flux   = energy_flux * scalar;

   return result;
}

// Numerical flux function
void Material::num_flux (const PrimVar& left,
                         const PrimVar& right,
                         const Vector& normal,
                         Flux& flux)
{
   
   double area = normal.norm();
   Vector unit_normal = normal / area;

   // Enthalpy
   double h_left  = GAMMA*left.pressure/(left.density*(GAMMA-1.0))
      + 0.5 * left.velocity.square();
   double h_right = GAMMA*right.pressure/(right.density*(GAMMA-1.0))
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
   double c = sqrt( (GAMMA-1.0) * (h - 0.5*velocity.square()) );

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

// Flux on slip walls
Flux Material::slip_flux (const PrimVar& state,
                          const Vector& normal)
{
   Flux flux;

   flux.mass_flux     = 0.0;
   flux.momentum_flux = normal * state.pressure;
   flux.energy_flux   = 0.0;

   return flux;
}
