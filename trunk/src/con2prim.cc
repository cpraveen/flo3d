#include "material.h"

PrimVar Material::con2prim (const ConVar& con_var)
{
   PrimVar prim_var;

   prim_var.density  = con_var.density;
   prim_var.velocity = con_var.momentum / con_var.density;
   prim_var.pressure = (gamma - 1.0) * 
        ( con_var.energy - 0.5 * con_var.momentum.square() / con_var.density );

   return prim_var;
}
