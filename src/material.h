#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include <string>
#include <cmath>
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
      PrimVar  operator*  (const PrimVar& prim_var) const; // componentwise multi
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
      double prandtl;
      double Cp;
      double T_0, T_ref, mu_ref; // constants for sutherland law
      std::string model;
      std::string flux;

      void initialize ();
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
                          const Vector&  normal,
                          Flux& flux);
      void viscous_flux (const PrimVar& state, 
                         const Vector&  dU, 
                         const Vector&  dV, 
                         const Vector&  dW, 
                         const Vector&  dT, 
                         const Vector&  normal, 
                         Flux&          flux);
      double viscosity (const double T);
      double temperature (const PrimVar& state);

};

inline
ConVar& ConVar::operator= (const ConVar& con_var)
{
   density  = con_var.density;
   momentum = con_var.momentum;
   energy   = con_var.energy;

   return *this;
}

inline
ConVar& ConVar::operator+= (const ConVar& con_var)
{
   density  += con_var.density;
   momentum += con_var.momentum;
   energy   += con_var.energy;

   return *this;
}

inline
ConVar ConVar::operator+ (const ConVar& con_var) const
{
   ConVar result;

   result.density  = density  + con_var.density;
   result.momentum = momentum + con_var.momentum;
   result.energy   = energy   + con_var.energy;

   return result;
}

inline
ConVar ConVar::operator- (const ConVar& con_var) const
{
   ConVar result;

   result.density  = density  - con_var.density;
   result.momentum = momentum - con_var.momentum;
   result.energy   = energy   - con_var.energy;

   return result;
}

inline
ConVar ConVar::operator- (const Flux& flux) const
{
   ConVar result;

   result.density  = density  - flux.mass_flux;
   result.momentum = momentum - flux.momentum_flux;
   result.energy   = energy   - flux.energy_flux;

   return result;
}

inline
ConVar ConVar::operator* (const double scalar) const
{
   ConVar result;

   result.density  = density  * scalar;
   result.momentum = momentum * scalar; 
   result.energy   = energy   * scalar;

   return result;
}

// Add two primitive variables
inline
PrimVar PrimVar::operator+ (const PrimVar& prim_var) const
{
   PrimVar result;

   result.density  = density  + prim_var.density;
   result.velocity = velocity + prim_var.velocity;
   result.pressure = pressure + prim_var.pressure;

   return result;
}

// Subtract two primitive variables
inline
PrimVar PrimVar::operator- (const PrimVar& prim_var) const
{
   PrimVar result;

   result.density  = density  - prim_var.density;
   result.velocity = velocity - prim_var.velocity;
   result.pressure = pressure - prim_var.pressure;

   return result;
}

// Multiply primitive by scalar and return result
inline
PrimVar PrimVar::operator* (const double& scalar) const
{
   PrimVar result;

   result.density  = density  * scalar;
   result.velocity = velocity * scalar; 
   result.pressure = pressure * scalar;

   return result;
}

// Multiply two primitive variables componentwise
// Result is another primitive variable
// NOTE: This is not scalar dot product
inline
PrimVar PrimVar::operator* (const PrimVar& prim_var) const
{
   PrimVar result;

   result.density  = density  * prim_var.density;
   result.velocity = velocity * prim_var.velocity;
   result.pressure = pressure * prim_var.pressure;

   return result;
}

// Multiply given primitive by scalar
inline
PrimVar& PrimVar::operator*= (const double& scalar)
{
   density  *= scalar;
   velocity *= scalar; 
   pressure *= scalar;

   return *this;
}

// Add primitive variable to given primitive variable
inline
PrimVar& PrimVar::operator+= (const PrimVar& prim_var)
{
   density  += prim_var.density;
   velocity += prim_var.velocity;
   pressure += prim_var.pressure;

   return *this;
}
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
double Material::viscosity (const double T)
{
   return mu_ref * std::pow(T/T_ref, 1.5) * (T_ref + T_0) / (T + T_0);
}

//------------------------------------------------------------------------------
//  Compute temperature given primitive state
//------------------------------------------------------------------------------
inline
double Material::temperature (const PrimVar& state)
{
   return state.pressure / (gas_const * state.density);
}

#endif
