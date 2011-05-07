#ifndef __FPARSER_H__
#define __FPARSER_H__

#include "fparser.hh"
#include <cmath>

class FParser: public FunctionParser
{
   public:
      FParser()
      {
         AddConstant("pi", M_PI);
      }
};

#endif
