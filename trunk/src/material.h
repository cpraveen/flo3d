#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include "vec.h"

const double GAMMA = 1.4;

// Primitive variable
class PrimVar
{
   public:
      double density, pressure;
      Vec    velocity;
};

class Flux
{
   public:
      double mass_flux;
      Vec    momentum_flux;
      double energy_flux;

      Flux& operator+= (const Flux& flux);
      Flux& operator-= (const Flux& flux);
      Flux& operator*= (const double& scalar);
      Flux  operator+  (const Flux& flux);
      Flux  operator*  (const double scalar);

      void zero ();
};

// Conserved variable
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
      Vec    momentum;

      void zero ();
      double pressure () const;

};

class Material
{
   public:

      ConVar  prim2con (const PrimVar& prim_var);
      PrimVar con2prim (const ConVar&  con_var);
      void    num_flux (const ConVar&, const ConVar&, const Vec&, Flux&);
      Flux    slip_flux (const ConVar& state, const Vec& normal);
};

#endif
