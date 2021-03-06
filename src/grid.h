#ifndef __GRID_H__
#define __GRID_H__

#include <vector>
#include "parameter.h"
#include "vec.h"
#include "face.h"

class Cell
{
   public:
      Vector       centroid;
      unsigned int vertex[4];
      int          face[4];
      double       volume;
      double       weight[4];
};

class Grid
{
   public:
      Grid () { n_vertex = n_cell = n_face = n_boundary_face = 0; };
      unsigned int n_vertex;
      unsigned int n_cell;
      unsigned int n_face;
      unsigned int n_boundary_face;
      double min_cell_volume;
      double max_cell_volume;
      double h_min, h_max;
      std::vector<Vector> vertex;
      std::vector<Cell>   cell;
      std::vector<Face>   face;
      std::vector<double> vertex_volume;

      void read (const Parameter& param);
      void find_cell_neighbour(const unsigned int& face_no,
                               const unsigned int& cell_no,
                               int&                neighbour_cell_no);


   private:
      void read_gmsh (std::string grid_file);
      void check_face_type (const std::map<int,BoundaryCondition>& bc);
      void preproc ();
      void compute_cell_centroid ();
      void compute_face_centroid ();
      void compute_cell_volume ();
      void compute_face_normal_and_area ();
      void compute_h ();
      void add_face (const Face& new_face);
      void make_faces ();
      void weight_average () ;
      void vertex_weight_check () ;
      void find_cell_faces ();
      void info ();
      void renumber_cell();

      std::vector< std::vector<unsigned int> > node_face;

};

#endif
