#include <cmath>
#include "material.h"

using namespace std;

//------------------------------------------------------------------------------
// KFVS split flux
//------------------------------------------------------------------------------
void Material::kfvs_split_flux (const double   sign,
                                const Vector&  normal,
                                const PrimVar& state,
                                Flux&          flux)
{
   // Normal velocity
   const double un = state.velocity * normal;

   const double beta = 0.5 * state.density / state.pressure;
   const double s    = un * sqrt(beta);
   const double A    = 0.5 * (1.0 + sign * erf(s));
   const double B    = 0.5 * sign * exp(-s*s) / sqrt(M_PI * beta);
   const double uf   = un * A + B;
   const double E    = total_energy (state);

   flux.mass_flux = state.density * uf;
   flux.momentum_flux = normal * state.pressure * A + 
                        state.velocity * state.density * uf;
   flux.energy_flux = (E + state.pressure) * A * un +
                      (E + 0.5 * state.pressure) * B;
}

//------------------------------------------------------------------------------
// KFVS flux function
//------------------------------------------------------------------------------
void Material::kfvs (const PrimVar& left,
                     const PrimVar& right,
                     const Vector& normal,
                     Flux& flux)
{
   Flux pos, neg;

   double area = normal.norm();
   Vector unit_normal = normal / area;

   kfvs_split_flux (+1, unit_normal, left,  pos);
   kfvs_split_flux (-1, unit_normal, right, neg);
   
   flux  = pos + neg;
   flux *= area;
}
