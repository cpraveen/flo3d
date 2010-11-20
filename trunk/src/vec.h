#ifndef __VEC_H__
#define __VEC_H__

// 3-D vector
class Vec
{
   public:
      double x, y, z;
      Vec&    operator=  (const Vec rhs);
      Vec&    operator=  (const double scalar);
      Vec&    operator+= (const Vec rhs);
      Vec&    operator-= (const Vec rhs);
      Vec     operator/  (const double scalar) const;
      Vec     operator*  (const double scalar) const;
      double  operator*  (const Vec vec) const; // Dot product of two vectors
      Vec     operator+  (const Vec vec) const;
      Vec     operator-  (const Vec vec) const;
      Vec     operator^  (const Vec vec) const; // Cross product of two vectors
      double  square () const;
      double  norm () const;
};

#endif
