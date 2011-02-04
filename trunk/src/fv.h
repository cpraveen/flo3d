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

      std::vector<PrimVar> primitive;
      std::vector<ConVar>  conserved_old;
      std::vector<PrimVar> primitive_vertex;
      std::vector<Flux>    residual;
      Flux                 residual_norm;
      double               residual_norm_total;
      double               residual_norm_total0;
      std::vector<double>  dt;
      double               dt_global;

      void reconstruct (const unsigned int,
                        const unsigned int,
                        const unsigned int,
                        const unsigned int,
                        std::vector<PrimVar>&) const;

      void reconstruct (const unsigned int vl,
                        const unsigned int cl,
                        PrimVar&           state) const;

      void initialize ();
      void interpolate_vertex ();
      void store_conserved_old ();
      void compute_residual ();
      void compute_dt ();
      void compute_residual_norm (const unsigned int iter);
      void update_solution (const unsigned int r);
      void solve ();
      void output (const unsigned int iter);
};

#endif
