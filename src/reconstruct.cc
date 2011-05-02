#include <cmath>
#include "fv.h"

#define LIMITER(r)  (2*(r)/(1+(r)*(r)))

using namespace std;

//------------------------------------------------------------------------------
// First order Reconstruct left and right states
//------------------------------------------------------------------------------
void FiniteVolume::reconstruct_first
(
 const unsigned int& f,
 bool                has_right,
 vector<PrimVar>&    state
) const
{
   // Left state
   unsigned int cl = grid.face[f].lcell;
   state[0] = primitive[cl];

   // Right state
   if(has_right)
   {
      unsigned int cr = grid.face[f].rcell;
      state[1] = primitive[cr];
   }
}

//------------------------------------------------------------------------------
// Second order Reconstruct left and right states
//------------------------------------------------------------------------------
void FiniteVolume::reconstruct_second
(
 const unsigned int& f,
 bool                has_right,
 vector<PrimVar>&    state
) const
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
   state[0] = primitive[cl] + ( face_avg - primitive_vertex[vl] ) * 0.25;

   // Right state
   if(has_right)
   {
      unsigned int vr = grid.face[f].rvertex;
      unsigned int cr = grid.face[f].rcell;
      state[1] = primitive[cr] + ( face_avg - primitive_vertex[vr] ) * 0.25;
   }
}

//------------------------------------------------------------------------------
// Apply limiter
//------------------------------------------------------------------------------
void limited_state(const PrimVar& prim_v, /* vertex value */
                   const PrimVar& prim_c, /* cell value */
                   const PrimVar& prim_f, /* face value opposite to vertex */
                   const PrimVar& prim_n, /* neighbour cell value */
                         PrimVar& state)  /* reconstructed state */
{
      const PrimVar dPrim = (prim_n - prim_c) * 0.5;
      const PrimVar dCV   = (prim_c - prim_v) * (1.0/3.0);
      const PrimVar fact  = dCV * dPrim;

      // First order reconstruction
      state = prim_c;

      // Now we add corrections with limiting

      // density
      if(fact.density > 0.0)
      {
         double r   = dCV.density / dPrim.density;
         double phi = LIMITER (r);
         state.density += phi * dCV.density;
      }

      // x velocity
      if(fact.velocity.x > 0.0)
      {
         double r   = dCV.velocity.x / dPrim.velocity.x;
         double phi = LIMITER (r);
         state.velocity.x += phi * dCV.velocity.x;
      }

      // y velocity
      if(fact.velocity.y > 0.0)
      {
         double r   = dCV.velocity.y / dPrim.velocity.y;
         double phi = LIMITER (r);
         state.velocity.y += phi * dCV.velocity.y;
      }

      // z velocity
      if(fact.velocity.z > 0.0)
      {
         double r   = dCV.velocity.z / dPrim.velocity.z;
         double phi = LIMITER (r);
         state.velocity.z += phi * dCV.velocity.z;
      }

      // pressure
      if(fact.pressure > 0.0)
      {
         double r   = dCV.pressure / dPrim.pressure;
         double phi = LIMITER (r);
         state.pressure += phi * dCV.pressure;
      }
}

//------------------------------------------------------------------------------
// Reconstruct left and right states
//------------------------------------------------------------------------------
void FiniteVolume::reconstruct_limited
(
 const unsigned int& f,
 bool                has_right,
 vector<PrimVar>&    state
) const
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

      PrimVar dCV = primitive[cl] - primitive_vertex[vl];

      // First order reconstruction
      state[0] = primitive[cl] + dCV * (1.0/3.0);
   }
}
