#include <iostream>
#include <algorithm>
#include <cmath>
#include "material.h"

using namespace std;

// Set conserved variable to zero
void PrimVar::zero ()
{
   density  = 0.0;
   velocity = 0.0;
   pressure = 0.0;
}

//------------------------------------------------------------------------------
// Do some initializations
//------------------------------------------------------------------------------
void Material::initialize ()
{
   Cp = gamma * gas_const / (gamma - 1.0);
}

//------------------------------------------------------------------------------
// Numerical flux function
//------------------------------------------------------------------------------
void Material::num_flux (const PrimVar& left,
                         const PrimVar& right,
                         const Vector& normal,
                         Flux& flux) const
{
   switch (flux_scheme)
   {
      case roe:
         roe_flux (left, right, normal, flux);
         break;

      case kfvs:
         kfvs_flux (left, right, normal, flux);
         break;

      default:
         cout << "num_flux: unknown flux " << flux_scheme << endl;
         abort ();
   }
}

//------------------------------------------------------------------------------
// Flux on slip walls
//------------------------------------------------------------------------------
void Material::slip_flux (const PrimVar& state,
                          const Vector&  normal,
                          Flux&          flux) const
{
   flux.mass_flux     = 0.0;
   flux.momentum_flux = normal * state.pressure;
   flux.energy_flux   = 0.0;
}

//------------------------------------------------------------------------------
// Euler Flux Calculation
//------------------------------------------------------------------------------
void Material::euler_flux (const PrimVar& prim, 
                           const Vector&  normal,
                           Flux&          flux) const
{
   // Enthalpy 
   double h  = gamma * prim.pressure / (prim.density * (gamma-1.0)) + 
               0.5 * prim.velocity.square();
   // Normal velocity
   double vn = prim.velocity * normal;

   flux.mass_flux = prim.density * vn;
   flux.momentum_flux = normal * prim.pressure + 
                        prim.velocity * flux.mass_flux;
   flux.energy_flux = h * flux.mass_flux;
}

//------------------------------------------------------------------------------
// viscous flux: TODO check formulae
//------------------------------------------------------------------------------
void Material::viscous_flux (const PrimVar& state, 
                             const Vector&  dU,
                             const Vector&  dV,
                             const Vector&  dW,
                             const Vector&  dT,
                             const Vector&  normal,
                             Flux&          flux) const
{
   double T = temperature (state);
   double mu = viscosity (T);
   double k = mu * Cp / prandtl;

   // Heat flux
   double q = -k * (dT * normal);

   // Divergence of velocity
   double div = dU.x + dV.y + dW.z;

   // Shear stress tensor: symmetric, compute only upper part
   double sxx = 2.0 * mu * (dU.x - (1.0/3.0) * div);
   double syy = 2.0 * mu * (dV.y - (1.0/3.0) * div);
   double szz = 2.0 * mu * (dW.z - (1.0/3.0) * div);
   double sxy = mu * (dU.y + dV.x);
   double sxz = mu * (dU.z + dW.x);
   double syz = mu * (dV.z + dW.y);

   flux.mass_flux = 0.0;
   flux.momentum_flux.x = -(sxx * normal.x + sxy * normal.y + sxz * normal.z);
   flux.momentum_flux.y = -(sxy * normal.x + syy * normal.y + syz * normal.z);
   flux.momentum_flux.z = -(sxz * normal.x + syz * normal.y + szz * normal.z);
   flux.energy_flux = flux.momentum_flux * state.velocity + q;
}
