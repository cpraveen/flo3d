#ifndef __FACE_H__
#define __FACE_H__

#include "vec.h"

class Face
{
   public:
      unsigned int vertex[3];
      int          lcell, rcell;
      int          lvertex, rvertex;
      Vector       normal;
      Vector       centroid;
      int          type;
      double       area;

      bool operator== (const Face& face) const;
};

inline
bool Face::operator== (const Face& face) const
{

   if ( vertex[0]==face.vertex[0] &&
        vertex[1]==face.vertex[1] &&
        vertex[2]==face.vertex[2]) return true;

   if ( vertex[0]==face.vertex[0] &&
        vertex[1]==face.vertex[2] &&
        vertex[2]==face.vertex[1]) return true;

   if ( vertex[0]==face.vertex[1] &&
        vertex[1]==face.vertex[2] &&
        vertex[2]==face.vertex[0]) return true;

   if ( vertex[0]==face.vertex[1] &&
        vertex[1]==face.vertex[0] &&
        vertex[2]==face.vertex[2]) return true;

   if ( vertex[0]==face.vertex[2] &&
        vertex[1]==face.vertex[0] &&
        vertex[2]==face.vertex[1]) return true;

   if ( vertex[0]==face.vertex[2] &&
        vertex[1]==face.vertex[1] &&
        vertex[2]==face.vertex[0]) return true;

   return false;
}

#endif
