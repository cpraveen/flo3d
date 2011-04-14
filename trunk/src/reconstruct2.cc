#include <cmath>
#include "fv.h"

#define LIMITER(r)  (2*(r)/(1+(r)*(r)))

using namespace std;

void limited_state(const PrimVar& prim_v, /* vertex value */
                   const PrimVar& prim_c, /* cell value */
                   const PrimVar& prim_f, /* face value opposite to vertex */
                   const PrimVar& prim_n, /* neighbour cell value */
                         PrimVar& state)  /* reconstructed state */
{
      PrimVar  dPrim = prim_n - prim_c;

      // left state reconstruction
      PrimVar dFV = prim_f - prim_v;
      PrimVar dCV = prim_c - prim_v;

      PrimVar fact1 = dFV * dCV;
      PrimVar fact2 = dFV * dPrim;

      // left density
      if(fact1.density > 0.0 && fact2.density > 0.0)
      {
         double r   = 0.5 * dFV.density / dPrim.density;
         double phi = LIMITER (r);
         state.density = prim_c.density + 
                         phi * 0.25 * dFV.density;
      }
      else
         state.density = prim_c.density;

      // left x velocity
      if(fact1.velocity.x > 0.0 && fact2.velocity.x > 0.0)
      {
         double r   = 0.5 * dFV.velocity.x / dPrim.velocity.x;
         double phi = LIMITER (r);
         state.velocity.x = prim_c.velocity.x + 
                            phi * 0.25 * dFV.velocity.x;
      }
      else
         state.velocity.x = prim_c.velocity.x;

      // left y velocity
      if(fact1.velocity.y > 0.0 && fact2.velocity.y > 0.0)
      {
         double r   = 0.5 * dFV.velocity.y / dPrim.velocity.y;
         double phi = LIMITER (r);
         state.velocity.y = prim_c.velocity.y + 
                            phi * 0.25 * dFV.velocity.y;
      }
      else
         state.velocity.y = prim_c.velocity.y;

      // left z velocity
      if(fact1.velocity.z > 0.0 && fact2.velocity.z > 0.0)
      {
         double r   = 0.5 * dFV.velocity.z / dPrim.velocity.z;
         double phi = LIMITER (r);
         state.velocity.z = prim_c.velocity.z + 
                               phi * 0.25 * dFV.velocity.z;
      }
      else
         state.velocity.z = prim_c.velocity.z;

      // left pressure
      if(fact1.pressure > 0.0 && fact2.pressure > 0.0)
      {
         double r   = 0.5 * dFV.pressure / dPrim.pressure;
         double phi = LIMITER (r);
         state.pressure = prim_c.pressure + 
                          phi * 0.25 * dFV.pressure;
      }
      else
         state.pressure = prim_c.pressure;
}

//------------------------------------------------------------------------------
// Reconstruct left and right states
//------------------------------------------------------------------------------
void FiniteVolume::reconstruct2(const unsigned int& f,
                                bool                has_right,
                                vector<PrimVar>&    state) const
{
   // Average on face
   PrimVar face_avg;
   face_avg.zero ();

   for(unsigned int i=0; i<3; ++i)
      face_avg += primitive_vertex[ grid.face[f].vertex[i] ];
   face_avg *= (1.0/3.0);

   // Left state
   unsigned int vl = grid.face[f].lvertex;
   unsigned int cl = grid.face[f].lcell;

   // For interior faces, get both states
   if(has_right)
   {
      unsigned int vr = grid.face[f].rvertex;
      unsigned int cr = grid.face[f].rcell;

      // left state
      limited_state (primitive_vertex[vl],
                     primitive[cl],
                     face_avg,
                     primitive[cr],
                     state[0]);

      // right state
      limited_state (primitive_vertex[vr],
                     primitive[cr],
                     face_avg,
                     primitive[cl],
                     state[1]);

   }
   else
   {
      // Boundary face
      // Only left state is present

      PrimVar dFV = face_avg - primitive_vertex[vl];
      PrimVar dCV = primitive[cl] - primitive_vertex[vl];
      PrimVar fact = dFV * dCV;

      // density
      if(fact.density > 0.0)
         state[0].density = primitive[cl].density + 
            0.25 * dFV.density;
      else
         state[0].density = primitive[cl].density;

      // x velocity
      if(fact.velocity.x > 0.0)
         state[0].velocity.x = primitive[cl].velocity.x + 
            0.25 * dFV.velocity.x;
      else
         state[0].velocity.x = primitive[cl].velocity.x;

      // y velocity
      if(fact.velocity.y > 0.0)
         state[0].velocity.y = primitive[cl].velocity.y + 
            0.25 * dFV.velocity.y;
      else
         state[0].velocity.y = primitive[cl].velocity.y;

      // z velocity
      if(fact.velocity.z > 0.0)
         state[0].velocity.z = primitive[cl].velocity.z + 
            0.25 * dFV.velocity.z;
      else
         state[0].velocity.z = primitive[cl].velocity.z;

      // pressure
      if(fact.pressure > 0.0)
         state[0].pressure = primitive[cl].pressure + 
            0.25 * dFV.pressure;
      else
         state[0].pressure = primitive[cl].pressure;

   }
}
