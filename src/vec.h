#ifndef __VEC_H__
#define __VEC_H__

// 3-D vector
class Vector
{
   public:
      double x, y, z;
      Vector& operator=  (const Vector rhs);
      Vector& operator=  (const double scalar);
      Vector& operator+= (const Vector rhs);
      Vector& operator-= (const Vector rhs);
      Vector& operator*= (const double scalar);
      Vector  operator/  (const double scalar) const;
      Vector  operator*  (const double scalar) const;
      double  operator*  (const Vector vec) const; // Dot product of two vectors
      Vector  operator+  (const Vector vec) const;
      Vector  operator-  (const Vector vec) const;
      Vector  operator^  (const Vector vec) const; // Cross product of two vectors
      double  square () const;
      double  norm () const;
};

#endif
