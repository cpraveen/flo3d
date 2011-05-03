#ifndef __FLUX_H__
#define __FLUX_H__

#include "vec.h"

//------------------------------------------------------------------------------
// Flux variable
//------------------------------------------------------------------------------
class Flux
{
   public:
      double mass_flux;
      Vector momentum_flux;
      double energy_flux;

      Flux& operator+= (const Flux& flux);
      Flux& operator-= (const Flux& flux);
      Flux& operator*= (const double& scalar);
      Flux  operator+  (const Flux& flux);
      Flux  operator*  (const double scalar);
      Flux  operator-  (const Flux& flux);
      void zero ();
};

//------------------------------------------------------------------------------
// Set all flux components to zero
//------------------------------------------------------------------------------
inline
void Flux::zero ()
{
   mass_flux     = 0.0;
   momentum_flux = 0.0;
   energy_flux   = 0.0;
}

//------------------------------------------------------------------------------
// Add flux to given flux
//------------------------------------------------------------------------------
inline
Flux& Flux::operator+= (const Flux& flux)
{
   mass_flux     += flux.mass_flux;
   momentum_flux += flux.momentum_flux;
   energy_flux   += flux.energy_flux;

   return *this;
}

//------------------------------------------------------------------------------
// Subtract flux from given flux
//------------------------------------------------------------------------------
inline
Flux& Flux::operator-= (const Flux& flux)
{
   mass_flux     -= flux.mass_flux;
   momentum_flux -= flux.momentum_flux;
   energy_flux   -= flux.energy_flux;

   return *this;
}

//------------------------------------------------------------------------------
// Multiply given flux by a scalar
//------------------------------------------------------------------------------
inline
Flux& Flux::operator*= (const double& scalar)
{
   mass_flux     *= scalar;
   momentum_flux *= scalar;
   energy_flux   *= scalar;

   return *this;
}

//------------------------------------------------------------------------------
// Add two fluxes
//------------------------------------------------------------------------------
inline
Flux Flux::operator+ (const Flux& flux)
{
   Flux result;

   result.mass_flux     = mass_flux + flux.mass_flux;
   result.momentum_flux = momentum_flux + flux.momentum_flux;
   result.energy_flux   = energy_flux + flux.energy_flux;

   return result;
}

//------------------------------------------------------------------------------
// Subtract two fluxes
//------------------------------------------------------------------------------
inline
Flux Flux::operator- (const Flux& flux)
{
   Flux result;
   result.mass_flux     = mass_flux - flux.mass_flux;
   result.momentum_flux = momentum_flux - flux.momentum_flux;
   result.energy_flux   = energy_flux - flux.energy_flux;
   return result;
}   

//------------------------------------------------------------------------------
// Multiply flux by a scalar and return result
//------------------------------------------------------------------------------
inline
Flux Flux::operator* (const double scalar)
{
   Flux result;

   result.mass_flux     = mass_flux * scalar;
   result.momentum_flux = momentum_flux * scalar;
   result.energy_flux   = energy_flux * scalar;

   return result;
}

#endif
