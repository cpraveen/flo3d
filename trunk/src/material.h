#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include <string>
#include "vec.h"

const double GAMMA = 1.4;

// Primitive variable
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
      Vector momentum;

      double pressure () const;

};

class Material
{
   public:
      double gamma;
      double gas_const;
      std::string model;
      std::string flux;

      ConVar  prim2con (const PrimVar& prim_var);
      PrimVar con2prim (const ConVar&  con_var);
      void    num_flux (const PrimVar& left, const PrimVar& right, const Vector& normal, Flux& flux);
      void    slip_flux (const PrimVar& state, const Vector& normal, Flux& flux);
};

#endif
