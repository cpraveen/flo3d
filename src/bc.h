#ifndef __BC_H__
#define __BC_H__

#include <iostream>
#include <string>
#include <cassert>
#include "face.h"
#include "fparser.h"
#include "primvar.h"
#include "material.h"

namespace BC
{
   enum BCType { none, slip, noslip, farfield, inlet, outlet, pressure };
}

//------------------------------------------------------------------------------
// Boundary condition class
//------------------------------------------------------------------------------
class BoundaryCondition
{
   public:
      BoundaryCondition () {};
      BoundaryCondition (Material                 &material,
                         std::string              &bc_type,
                         std::vector<std::string> &variable,
                         std::vector<std::string> &function);
      void apply (const Face            &face,
                  std::vector<PrimVar> &state);
      void apply (const Vector  &vertex,
                  PrimVar       &state);
      void apply_slip (const Face           &face,
                       std::vector<PrimVar> &state);
      void apply_noslip (const Face           &face,
                         std::vector<PrimVar> &state);
      void apply_noslip (const Vector &vertex,
                         PrimVar      &state);
      void apply_pressure (const Face           &face,
                           std::vector<PrimVar> &state);
      void apply_pressure (const Vector &vertex,
                           PrimVar      &state);
      void apply_inlet (const Face           &face,
                        std::vector<PrimVar> &state);
      void apply_inlet (const Vector  &vertex,
                        PrimVar       &state);
      void apply_outlet (const Face           &face,
                         std::vector<PrimVar> &state);
      void apply_farfield (const Face           &face,
                           std::vector<PrimVar> &state);
      std::string    name;
      BC::BCType     type;
      bool           adiabatic;

   private:
      Material*      material;
      FParser        density;
      FParser        xvelocity;
      FParser        yvelocity;
      FParser        zvelocity;
      FParser        pressure;
      FParser        temperature;
};

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
inline
BoundaryCondition::BoundaryCondition (Material                 &material,
                                      std::string              &bc_type,
                                      std::vector<std::string> &variable,
                                      std::vector<std::string> &function)
:
   name (bc_type),
   material(&material)
{
   // Set to none for safety purpose
   type = BC::none;

   // Slip bc, no state is required
   if(bc_type == "slip")
   {
      assert (variable.size() == 0);
      type = BC::slip;
   }
   // noslip bc: velocity is specified. If temperature is also specified
   // then it is also used. In this case, adiabatic bc is not used
   else if(bc_type == "noslip")
   {
      assert (variable.size() == 3 || variable.size() == 4);
      type = BC::noslip;
      adiabatic = true;
      bool has_xvelocity   = false;
      bool has_yvelocity   = false;
      bool has_zvelocity   = false;
      for(unsigned int i=0; i<variable.size(); ++i)
      {
         if(variable[i] == "xvelocity")
         {
            has_xvelocity = true;
            xvelocity.Parse (function[i], "x,y,z");
         }
         else if(variable[i] == "yvelocity")
         {
            has_yvelocity = true;
            yvelocity.Parse (function[i], "x,y,z");
         }
         else if(variable[i] == "zvelocity")
         {
            has_zvelocity = true;
            zvelocity.Parse (function[i], "x,y,z");
         }
         else if(variable[i] == "temperature")
         {
            temperature.Parse (function[i], "x,y,z");
            adiabatic = false;
         }
      }
      assert (has_xvelocity && has_yvelocity && has_zvelocity);
   }
   // In this case only pressure is specified
   else if(bc_type == "pressure")
   {
      assert (variable.size() == 1);
      type = BC::pressure;
      assert (variable[0] == "pressure");
      pressure.Parse (function[0], "x,y,z");
   }
   // All values are specified
   else if(bc_type == "inlet" || bc_type == "farfield")
   {
      assert (variable.size() == 5);
      if(bc_type == "inlet")
         type = BC::inlet;
      else
         type = BC::farfield;
      bool has_density   = false;
      bool has_xvelocity = false;
      bool has_yvelocity = false;
      bool has_zvelocity = false;
      bool has_pressure  = false;
      for(unsigned int i=0; i<variable.size(); ++i)
      {
         if(variable[i] == "density")
         {
            has_density = true;
            density.Parse (function[i], "x,y,z");
         }
         else if(variable[i] == "xvelocity")
         {
            has_xvelocity = true;
            xvelocity.Parse (function[i], "x,y,z");
         }
         else if(variable[i] == "yvelocity")
         {
            has_yvelocity = true;
            yvelocity.Parse (function[i], "x,y,z");
         }
         else if(variable[i] == "zvelocity")
         {
            has_zvelocity = true;
            zvelocity.Parse (function[i], "x,y,z");
         }
         else if(variable[i] == "pressure")
         {
            has_pressure = true;
            pressure.Parse (function[i], "x,y,z");
         }
      }
      assert (has_density && has_xvelocity && has_yvelocity && has_zvelocity &&
              has_pressure);
   }
   // At outflow nothing is specified
   else if(bc_type == "outlet")
   {
      assert (variable.size() == 0);
      type = BC::outlet;
   }
   else
   {
      std::cout << "BoundaryCondition: Unknown boundary condition " << bc_type << std::endl;
      abort();
   }

   if(type == BC::none)
   {
      std::cout << "BoundaryCondition: unknown bc for " << bc_type << std::endl;
      abort();
   }
}

//------------------------------------------------------------------------------
// Normal velocity is zero
// Find state[1] so that average state has zero normal velocity
// u[1] = u[0] - 2 (u[0] * n) n
// Then u[1] * n = - u[0] * n
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_slip(const Face           &face,
                                   std::vector<PrimVar> &state)
{
   state[1] = state[0];
   Vector unit_normal = face.normal / face.area;
   state[1].velocity -= unit_normal * (state[1].velocity * unit_normal) * 2.0;
}

//------------------------------------------------------------------------------
// Velocity is specified
// If temperature is specified, then pressure is computed using the temperature
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_noslip(const Face           &face,
                                     std::vector<PrimVar> &state)
{
   state[1].density  = state[0].density;

   double point[3]  = {face.centroid.x, face.centroid.y, face.centroid.z};
   state[1].velocity.x = 2.0 * xvelocity.Eval(point) - state[0].velocity.x;
   state[1].velocity.y = 2.0 * yvelocity.Eval(point) - state[0].velocity.y;
   state[1].velocity.z = 2.0 * zvelocity.Eval(point) - state[0].velocity.z;

   if(adiabatic)
      state[1].pressure = state[0].pressure;
   else
   {
      double T = temperature.Eval(point);
      state[1].pressure = material->gas_const * state[1].density * T;
   }
}

//------------------------------------------------------------------------------
// Set velocity for noslip point
// Used in vertex averaging 
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_noslip(const Vector &vertex,
                                     PrimVar      &state)
{
   double point[3]  = {vertex.x, vertex.y, vertex.z};
   state.velocity.x = xvelocity.Eval(point);
   state.velocity.y = yvelocity.Eval(point);
   state.velocity.z = zvelocity.Eval(point);

   if(!adiabatic)
   {
      double T = temperature.Eval(point);
      state.pressure = material->gas_const * state.density * T;
   }
}

//------------------------------------------------------------------------------
// Reset pressure value. Other states remain same
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_pressure (const Face           &face,
                                        std::vector<PrimVar> &state)
{
   double point[3] = {face.centroid.x, face.centroid.y, face.centroid.z};

   state[1] = state[0];
   state[1].pressure  = pressure.Eval(point);
}

//------------------------------------------------------------------------------
// Reset pressure value. Other states remain same
// Used in vertex averaging 
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_pressure (const Vector  &vertex,
                                        PrimVar       &state)
{
   double point[3]  = {vertex.x, vertex.y, vertex.z};
   state.pressure  = pressure.Eval(point);
}

//------------------------------------------------------------------------------
// At inlet all values are specified
// Both states are set to inlet values
// e.g., supersonic inlet
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_inlet (const Face           &face,
                                     std::vector<PrimVar> &state)
{
   double point[3]  = {face.centroid.x, face.centroid.y, face.centroid.z};
   state[1].density    = density.Eval(point);
   state[1].velocity.x = xvelocity.Eval(point);
   state[1].velocity.y = yvelocity.Eval(point);
   state[1].velocity.z = zvelocity.Eval(point);
   state[1].pressure   = pressure.Eval(point);

   state[0] = state[1];
}

//------------------------------------------------------------------------------
// At inlet all values are specified
// Used in vertex averaging 
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_inlet (const Vector &vertex,
                                     PrimVar      &state)
{
   double point[3]  = {vertex.x, vertex.y, vertex.z};
   state.density    = density.Eval(point);
   state.velocity.x = xvelocity.Eval(point);
   state.velocity.y = yvelocity.Eval(point);
   state.velocity.z = zvelocity.Eval(point);
   state.pressure   = pressure.Eval(point);
}

//------------------------------------------------------------------------------
// At outlet all all values are from inside values
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_outlet (const Face           &face,
                                      std::vector<PrimVar> &state)
{
   state[1] = state[0];
}

//------------------------------------------------------------------------------
// At farfield all values are specified
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_farfield (const Face           &face,
                                        std::vector<PrimVar> &state)
{
   double point[3]  = {face.centroid.x, face.centroid.y, face.centroid.z};
   state[1].density    = density.Eval(point);
   state[1].velocity.x = xvelocity.Eval(point);
   state[1].velocity.y = yvelocity.Eval(point);
   state[1].velocity.z = zvelocity.Eval(point);
   state[1].pressure   = pressure.Eval(point);
}

//------------------------------------------------------------------------------
// Apply boundary condition based on type
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply(const Face           &face,
                              std::vector<PrimVar> &state)
{
   switch(type)
   {
      case BC::slip:
         apply_slip (face, state);
         break;

      case BC::noslip:
         apply_noslip (face, state);
         break;

      case BC::pressure:
         apply_pressure (face, state);
         break;

      case BC::inlet:
         apply_inlet (face, state);
         break;

      case BC::outlet:
         apply_outlet (face, state);
         break;

      case BC::farfield:
         apply_farfield (face, state);
         break;

      default:
         std::cout << "BoundaryCondition::apply" << std::endl;
         std::cout << "   Unknown boundary condition: " << name << std::endl;
         abort ();
   }
}

//------------------------------------------------------------------------------
// Apply dirichlet boundary condition based on type
// This is for boundary vertices
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply(const Vector  &vertex,
                              PrimVar       &state)
{
   switch(type)
   {
      // Nothing to do
      case BC::slip:
         break;

      case BC::noslip:
         apply_noslip (vertex, state);
         break;

      case BC::pressure:
         apply_pressure (vertex, state);
         break;

      case BC::inlet:
         apply_inlet (vertex, state);
         break;

      // Nothing to do
      case BC::outlet:
         break;

      // Nothing to do
      case BC::farfield:
         break;

      default:
         std::cout << "BoundaryCondition::apply" << std::endl;
         std::cout << "   Unknown boundary condition: " << name << std::endl;
         abort ();
   }
}

#endif
