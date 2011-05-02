#ifndef __IC_H__
#define __IC_H__

#include <iostream>
#include <string>
#include "vec.h"
#include "material.h"
#include "fparser.hh"

//------------------------------------------------------------------------------
// Class to store initial condition functions
//------------------------------------------------------------------------------
class InitialCondition
{
   public:
      void    add (std::string, std::string);
      PrimVar value (const Vector& p);

   private:
      FunctionParser density;
      FunctionParser xvelocity;
      FunctionParser yvelocity;
      FunctionParser zvelocity;
      FunctionParser pressure;
};

//------------------------------------------------------------------------------
// Add function as defined in "fun"
//------------------------------------------------------------------------------
inline
void InitialCondition::add (std::string variable, std::string fun)
{
   if(variable == "density")
      density.Parse (fun, "x,y,z");
   else if(variable == "xvelocity")
      xvelocity.Parse (fun, "x,y,z");
   else if(variable == "yvelocity")
      yvelocity.Parse (fun, "x,y,z");
   else if(variable == "zvelocity")
      zvelocity.Parse (fun, "x,y,z");
   else if(variable == "pressure")
      pressure.Parse (fun, "x,y,z");
   else
   {
      std::cout << "InitialCondition::add: Unknown variable " << variable << std::endl;
      abort ();
   }
}

//------------------------------------------------------------------------------
// Evaluate primitive variables for given point
//------------------------------------------------------------------------------
inline
PrimVar InitialCondition::value (const Vector& p)
{
   PrimVar result;

   double vals[3] = {p.x, p.y, p.z};

   result.density    = density.Eval (vals);
   result.velocity.x = xvelocity.Eval (vals);
   result.velocity.y = yvelocity.Eval (vals);
   result.velocity.z = zvelocity.Eval (vals);
   result.pressure   = pressure.Eval (vals);

   return result;
}

#endif
