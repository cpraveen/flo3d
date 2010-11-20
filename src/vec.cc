#include <cmath>
#include "vec.h"

using namespace std;

// Assign one vector to another
Vec& Vec::operator= (const Vec rhs){
   x = rhs.x;
   y = rhs.y;
   z = rhs.z;

   return *this;
}

// Assign one vector to another
Vec& Vec::operator= (const double scalar){
   x = y = z = scalar;

   return *this;
}

// Assign one vector to another
Vec& Vec::operator+= (const Vec rhs){
   x += rhs.x;
   y += rhs.y;
   z += rhs.z;

   return *this;
}

// Add two vectors
Vec Vec::operator+  (const Vec vec) const
{
   Vec result;

   result.x = x + vec.x;
   result.y = y + vec.y;
   result.z = z + vec.z;

   return result;
}

// Subtract two vectors
Vec Vec::operator-  (const Vec vec) const
{
   Vec result;

   result.x = x - vec.x;
   result.y = y - vec.y;
   result.z = z - vec.z;

   return result;
}

// Divide a vector by a scalar
Vec Vec::operator/ (const double scalar) const
{
   Vec result;
   result.x = x / scalar;
   result.y = y / scalar;
   result.z = z / scalar;

   return result;
}

// Multiply a vector by a scalar
Vec Vec::operator* (const double scalar) const
{
   Vec result;
   result.x = x * scalar;
   result.y = y * scalar;
   result.z = z * scalar;

   return result;
}

// L2 norm square of vector
double Vec::square () const
{
   return x*x + y*y + z*z;
}

// L2 norm of vector
double Vec::norm () const
{
   return sqrt(x*x + y*y + z*z);
}

// Dot product of two vectors
double Vec::operator* (const Vec vec) const
{
   return x * vec.x + y * vec.y + z * vec.z;
}

// Cross product of two vectors
Vec Vec::operator^ (const Vec vec) const
{
   Vec result;

   result.x = y * vec.z - z * vec.y;
   result.y = z * vec.x - x * vec.z;
   result.z = x * vec.y - y * vec.x;

   return result;
}
