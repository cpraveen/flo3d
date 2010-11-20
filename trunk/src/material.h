#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include "vec.h"

// Primitive variable
class PrimVar
{
   public:
      double density, pressure;
      Vec    velocity;
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
      ConVar  operator*  (const double scalar) const;

      double density, energy;
      Vec    momentum;

      void zero ();

};

class Material
{
   public:
      static const unsigned int n_var = 5;
      static const double gamma = 1.4;

      ConVar  prim2con (const PrimVar& prim_var);
      PrimVar con2prim (const ConVar&  con_var);
};

#endif
