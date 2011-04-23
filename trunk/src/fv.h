#ifndef __FV_H__
#define __FV_H__

#include <fstream>
#include <vector>
#include "parameter.h"
#include "material.h"
#include "grid.h"

// Main class for finite volume scheme
class FiniteVolume
{
   public:
      FiniteVolume (char* file) 
      { 
         param.file = file;
         param.read ();
      }
      ~FiniteVolume () {};
      void run ();

   private:
      std::ofstream res;
      Parameter param;
      Grid      grid;

      std::vector<PrimVar> primitive;
      std::vector<ConVar>  conserved_old;
      std::vector<PrimVar> primitive_vertex;
      std::vector<Flux>    residual;
      std::vector<Vector>  dU, dV, dW, dT;
      Flux                 residual_norm;
      double               residual_norm_total;
      double               residual_norm_total0;
      std::vector<double>  dt;
      double               dt_global;
      double               elapsed_time;

      void reconstruct (const unsigned int&      f,
                        bool                     has_right,
                        std::vector<PrimVar>&    state) const;

      // Adding limited reconstruction
      void reconstruct2(const unsigned int&      f,
                        bool                     has_right,
                        std::vector<PrimVar>&    state) const;

      void initialize ();
      void interpolate_vertex ();
      void compute_vertex_gradients ();
      void store_conserved_old ();
      void compute_residual ();
      void compute_dt ();
      void compute_residual_norm (const unsigned int iter);
      void log_messages (const unsigned int iter);
      void update_solution (const unsigned int r);
      void solve ();
      void output (const unsigned int iter);
      void output_restart ();
      void lusgs ();
      
};

#endif
