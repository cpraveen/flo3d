#include "parameter.h"
#include "fv.h"

ConVar Material::prim2con(const PrimVar& prim_var)
{
   ConVar con_var;

   con_var.density  = prim_var.density;
   con_var.momentum = prim_var.velocity * prim_var.density;
   con_var.energy   = prim_var.pressure/(gamma - 1.0) +
                        0.5 * prim_var.velocity.square() * prim_var.density;

   return con_var;
}
