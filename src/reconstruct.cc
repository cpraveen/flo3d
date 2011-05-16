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
   // Left state
   unsigned int vl = grid.face[f].lvertex;
   unsigned int cl = grid.face[f].lcell;
   state[0] = primitive[cl] + ( primitive[cl] - primitive_vertex[vl] ) * (1.0/3.0);

   // Right state
   if(has_right)
   {
      unsigned int vr = grid.face[f].rvertex;
      unsigned int cr = grid.face[f].rcell;
      state[1] = primitive[cr] + ( primitive[cr] - primitive_vertex[vr] ) * (1.0/3.0);
   }
}

//------------------------------------------------------------------------------
// Second order Reconstruct left and right states
// Frink scheme
//------------------------------------------------------------------------------
void FiniteVolume::reconstruct_secondF
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
                     primitive[cr],
                     state[0]);

      // right state
      limited_state (primitive_vertex[vr],
                     primitive[cr],
                     primitive[cl],
                     state[1]);
   }
   else
   {
      // Boundary face
      // Only left state is present

      // Average on face
      PrimVar face_avg;
      face_avg.zero ();

      for(unsigned int i=0; i<3; ++i)
         face_avg += primitive_vertex[ grid.face[f].vertex[i] ];
      face_avg *= (1.0/3.0);

      PrimVar dCV = primitive[cl] - primitive_vertex[vl];
      PrimVar dFC = face_avg - primitive[cl];
      PrimVar fact= dCV * dFC;

      state[0] = primitive[cl];

      if(fact.density > 0.0)
         state[0].density += (1.0/3.0) * dCV.density;

      if(fact.velocity.x > 0.0)
         state[0].velocity.x += (1.0/3.0) * dCV.velocity.x;

      if(fact.velocity.y > 0.0)
         state[0].velocity.y += (1.0/3.0) * dCV.velocity.y;

      if(fact.velocity.z > 0.0)
         state[0].velocity.z += (1.0/3.0) * dCV.velocity.z;

      if(fact.pressure > 0.0)
         state[0].pressure += (1.0/3.0) * dCV.pressure;

   }
}

//------------------------------------------------------------------------------
// Apply limiter for Frink scheme
//------------------------------------------------------------------------------
void limitedF_state(const PrimVar& prim_v, /* vertex value */
                    const PrimVar& prim_c, /* cell value */
                    const PrimVar& prim_f, /* face value opposite to vertex */
                    const PrimVar& prim_n, /* neighbour cell value */
                          PrimVar& state)  /* reconstructed state */
{
      PrimVar  dPrim = prim_n - prim_c;

      // left state reconstruction
      PrimVar dFV = prim_f - prim_v;
      PrimVar dCV = prim_c - prim_v;

      PrimVar fact1 = dCV * dPrim;
      PrimVar fact2 = dFV * dPrim;

      // First order reconstruction
      state = prim_c;

      // Now we add corrections with limiting

      // left density
      if(fact1.density > 0.0 && fact2.density > 0.0)
      {
         double r   = 0.5 * dFV.density / dPrim.density;
         double phi = LIMITER (r);
         state.density += phi * 0.25 * dFV.density;
      }

      // left x velocity
      if(fact1.velocity.x > 0.0 && fact2.velocity.x > 0.0)
      {
         double r   = 0.5 * dFV.velocity.x / dPrim.velocity.x;
         double phi = LIMITER (r);
         state.velocity.x += phi * 0.25 * dFV.velocity.x;
      }

      // left y velocity
      if(fact1.velocity.y > 0.0 && fact2.velocity.y > 0.0)
      {
         double r   = 0.5 * dFV.velocity.y / dPrim.velocity.y;
         double phi = LIMITER (r);
         state.velocity.y += phi * 0.25 * dFV.velocity.y;
      }

      // left z velocity
      if(fact1.velocity.z > 0.0 && fact2.velocity.z > 0.0)
      {
         double r   = 0.5 * dFV.velocity.z / dPrim.velocity.z;
         double phi = LIMITER (r);
         state.velocity.z += phi * 0.25 * dFV.velocity.z;
      }

      // left pressure
      if(fact1.pressure > 0.0 && fact2.pressure > 0.0)
      {
         double r   = 0.5 * dFV.pressure / dPrim.pressure;
         double phi = LIMITER (r);
         state.pressure += phi * 0.25 * dFV.pressure;
      }
}

//------------------------------------------------------------------------------
// Reconstruct left and right states
// Frink scheme
//------------------------------------------------------------------------------
void FiniteVolume::reconstruct_limitedF
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
      limitedF_state (primitive_vertex[vl],
                      primitive[cl],
                      face_avg,
                      primitive[cr],
                      state[0]);

      // right state
      limitedF_state (primitive_vertex[vr],
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
      PrimVar dFC = face_avg - primitive[cl];
      PrimVar dCV = primitive[cl] - primitive_vertex[vl];
      PrimVar fact = dFC * dCV;

      // First order reconstruction
      state[0] = primitive[cl];

      // density
      if(fact.density > 0.0)
         state[0].density += 0.25 * dFV.density;

      // x velocity
      if(fact.velocity.x > 0.0)
         state[0].velocity.x += 0.25 * dFV.velocity.x;

      // y velocity
      if(fact.velocity.y > 0.0)
         state[0].velocity.y += 0.25 * dFV.velocity.y;

      // z velocity
      if(fact.velocity.z > 0.0)
         state[0].velocity.z += 0.25 * dFV.velocity.z;

      // pressure
      if(fact.pressure > 0.0)
         state[0].pressure += 0.25 * dFV.pressure;
   }
}

//------------------------------------------------------------------------------
// Reconstruct left and right states
//------------------------------------------------------------------------------
void FiniteVolume::reconstruct (const unsigned int& f,
                                bool                has_right,
                                vector<PrimVar>&    state) const
{
   switch(param.reconstruct_scheme)
   {
      case Parameter::first:
         reconstruct_first (f, has_right, state);
         break;

      case Parameter::second:
         reconstruct_second (f, has_right, state);
         break;

      case Parameter::secondF:
         reconstruct_secondF (f, has_right, state);
         break;

      case Parameter::limited:
         reconstruct_limited (f, has_right, state);
         break;

      case Parameter::limitedF:
         reconstruct_limitedF (f, has_right, state);
         break;

      default:
         cout << "reconstruct: unknown reconstruction scheme = " 
              << param.reconstruct_scheme << endl;
         abort ();
   }
}
