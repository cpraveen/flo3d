#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include <string>
#include "vec.h"

//------------------------------------------------------------------------------
// Primitive variable
//------------------------------------------------------------------------------
class PrimVar
{
   public:
      double density, pressure;
      Vector velocity;

      PrimVar  operator+  (const PrimVar& prim_var) const;
      PrimVar  operator-  (const PrimVar& prim_var) const;
      PrimVar  operator*  (const double& scalar) const;
      PrimVar& operator*= (const double& scalar);
      PrimVar& operator+= (const PrimVar& prim_var);
      void zero ();
};

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
// Conserved variable
//------------------------------------------------------------------------------
class ConVar
{
   public:
      ConVar () {};
      ~ConVar () {};
      ConVar& operator=  (const ConVar& con_var);
      ConVar& operator+= (const ConVar& con_var);
      ConVar  operator+  (const ConVar& con_var) const;
      ConVar  operator-  (const ConVar& con_var) const;
      ConVar  operator-  (const Flux&   flux) const;
      ConVar  operator*  (const double  scalar) const;

      double density, energy;
      Vector momentum;

};

//------------------------------------------------------------------------------
// Material class
//------------------------------------------------------------------------------
class Material
{
   public:
      double gamma;
      double gas_const;
      std::string model;
      std::string flux;

      ConVar  prim2con (const PrimVar& prim_var);
      PrimVar con2prim (const ConVar&  con_var);
      void    num_flux (const PrimVar& left, 
                        const PrimVar& right, 
                        const Vector& normal, 
                        Flux& flux);
      void    slip_flux (const PrimVar& state, 
                         const Vector& normal, 
                         Flux& flux);
      void    euler_flux (const PrimVar& prim, 
                          Flux& flux, 
                          const Vector& normal );

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

#endif
