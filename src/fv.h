#ifndef __FV_H__
#define __FV_H__

#include <vector>
#include "parameter.h"
#include "material.h"
#include "grid.h"

// Main class for finite volume scheme
class FiniteVolume
{
   public:
      FiniteVolume () {};
      ~FiniteVolume () {};
      void run ();

   private:
      Parameter param;
      Grid      grid;
      Material  material;

      std::vector<ConVar> solution;
      std::vector<ConVar> solution_old;
      std::vector<ConVar> solution_vertex;
      std::vector<ConVar> residual;
      std::vector<double> dt;

      void initialize ();
      void interpolate_vertex ();
      void compute_residual ();
      void update_solution ();
      void solve ();
      void output ();
};

#endif
