#ifndef __GRID_H__
#define __GRID_H__

#include <vector>
#include "parameter.h"
#include "vec.h"

class Cell
{
   public:
      unsigned int vertex[4];
      int          neighbour[4];
      double       volume;
      double       weight[4];
};

class Face
{
   public:
      unsigned int vertex[3];
      int          lcell, rcell;
      int          lvertex, rvertex;
      Vector       normal;
      int          type;

      bool operator== (const Face& face) const;
};

class Grid
{
   public:
      Grid () { n_vertex = n_cell = n_face = 0; };
      unsigned int n_vertex;
      unsigned int n_cell;
      unsigned int n_face;
      double min_cell_volume;
      double max_cell_volume;
      std::vector<Vector> vertex;
      std::vector<Cell>   cell;
      std::vector<Face>   face;

      void read (const Parameter& param);

   private:
      void read_gmsh (std::string grid_file);
      void check_face_type (const std::map<int,BCType>& bc_type);
      void preproc ();
      void compute_cell_volume ();
      void compute_face_normal ();
      void add_face (const Face& new_face);
      void make_faces ();
      void weight_average () ;
      void vertex_weight_check () ;
      void find_cell_surr_cell ();
      void info ();
      void renumbering_cell();

      std::vector< std::vector<unsigned int> > node_face;

};

#endif
