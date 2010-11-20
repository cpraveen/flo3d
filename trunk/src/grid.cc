#include "grid.h"

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
