#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include <iostream>
#include <string>
#include <cmath>
#include "vec.h"
#include "primvar.h"
#include "flux.h"
#include "convar.h"

//------------------------------------------------------------------------------
// Material class
//------------------------------------------------------------------------------
class Material
{
   public:
      double gamma;
      double gas_const;
      double prandtl;
      double Cp;
      double T_0, T_ref, mu_ref; // constants for sutherland law
      std::string model;
      enum FluxScheme { roe, kfvs };
      FluxScheme flux_scheme;

      enum MuModel {mu_constant, mu_sutherland};
      MuModel mu_model;

      void initialize ();
      ConVar  prim2con (const PrimVar& prim_var);
      PrimVar con2prim (const ConVar&  con_var);
      void num_flux(const PrimVar& left, 
                    const PrimVar& right, 
                    const Vector& normal, 
                    Flux& flux) const;
      void    roe_flux (const PrimVar& left, 
                        const PrimVar& right, 
                        const Vector& normal, 
                        Flux& flux) const;
      void kfvs_split_flux (const double   sign,
                            const Vector&  normal,
                            const PrimVar& state,
                            Flux&          flux) const;
      void    kfvs_flux(const PrimVar& left, 
                        const PrimVar& right, 
                        const Vector& normal, 
                        Flux& flux) const;
      void    slip_flux (const PrimVar& state, 
                         const Vector& normal, 
                         Flux& flux) const;
      void    euler_flux (const PrimVar& prim, 
                          const Vector&  normal,
                          Flux& flux) const;
      void viscous_flux (const PrimVar& state, 
                         const Vector&  dU, 
                         const Vector&  dV, 
                         const Vector&  dW, 
                         const Vector&  dT, 
                         const Vector&  normal, 
                         Flux&          flux) const;
      double viscosity (const double T) const;
      double temperature (const PrimVar& state) const;
      double total_energy (const PrimVar& state) const;
      double sound_speed (const PrimVar& state) const;

};

//------------------------------------------------------------------------------
// Convert primitive to conserved
//------------------------------------------------------------------------------
inline
ConVar Material::prim2con(const PrimVar& prim_var)
{
   ConVar con_var;

   con_var.density  = prim_var.density;
   con_var.momentum = prim_var.velocity * prim_var.density;
   con_var.energy   = prim_var.pressure/(gamma - 1.0) +
                        0.5 * prim_var.velocity.square() * prim_var.density;

   return con_var;
}

//------------------------------------------------------------------------------
// Convert conserved to primitive
//------------------------------------------------------------------------------
inline
PrimVar Material::con2prim (const ConVar& con_var)
{
   PrimVar prim_var;

   prim_var.density  = con_var.density;
   prim_var.velocity = con_var.momentum / con_var.density;
   prim_var.pressure = (gamma - 1.0) * 
        ( con_var.energy - 0.5 * con_var.momentum.square() / con_var.density );

   return prim_var;
}

//------------------------------------------------------------------------------
// Viscosity coefficient according to sutherland law
//------------------------------------------------------------------------------
inline
double Material::viscosity (const double T) const
{
   switch (mu_model)
   {
      case mu_constant:
         return mu_ref;

      case mu_sutherland:
         return mu_ref * std::pow(T/T_ref, 1.5) * (T_ref + T_0) / (T + T_0);

      default:
         std::cout << "viscosity: unknown model " << mu_model << std::endl;
         abort ();
   }
}

//------------------------------------------------------------------------------
//  Compute temperature given primitive state
//------------------------------------------------------------------------------
inline
double Material::temperature (const PrimVar& state) const
{
   return state.pressure / (gas_const * state.density);
}

//------------------------------------------------------------------------------
// Total energy per unit volume
//------------------------------------------------------------------------------
inline
double Material::total_energy (const PrimVar& state) const
{
   return state.pressure / (gamma - 1.0) + 
          0.5 * state.density * state.velocity.square();
}

//------------------------------------------------------------------------------
// sound speed
//------------------------------------------------------------------------------
inline
double Material::sound_speed (const PrimVar& state) const
{
   return sqrt(gamma * state.pressure / state.density);
}

#endif
