#include "material.h"

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

void ConVar::zero ()
{
   density  = 0.0;
   momentum = 0.0;
   energy   = 0.0;
}

double ConVar::pressure () const
{
   return (GAMMA - 1.0) * (energy - 0.5 * momentum.square() / density);
}

void Flux::zero ()
{
   mass_flux  = 0.0;
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

Flux Flux::operator* (const double scalar)
{
   Flux result;

   result.mass_flux     = mass_flux * scalar;
   result.momentum_flux = momentum_flux * scalar;
   result.energy_flux   = energy_flux * scalar;

   return result;
}

// Numerical flux function
Flux Material::num_flux (const ConVar& state_left,
                         const ConVar& state_right,
                         const Vec normal)
{
   Flux flux;

   return flux;
}

// Flux on slip walls
Flux Material::slip_flux (const ConVar& state,
                          const Vec& normal)
{
   Flux flux;

   double pressure = state.pressure ();

   flux.mass_flux     = 0.0;
   flux.momentum_flux = normal * pressure;
   flux.energy_flux   = 0.0;

   return flux;
}
