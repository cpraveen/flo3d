#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// For each force, create list of faces
//------------------------------------------------------------------------------
void FiniteVolume::create_force_face_list ()
{
   if(param.force_data.size() == 0)
      return;

   force.resize (param.force_data.size());

   // For each force, count how many faces of a given type were found
   // This is meant to catch a mistake where a user specifies a face type
   // which does not exist in the grid.
   vector< map<int,int> > nface(force.size());
   for(unsigned int j=0; j<force.size(); ++j)
   {
      ForceData& force_data = param.force_data[j];
      for(unsigned int k=0; k<force_data.face_type.size(); ++k)
         nface[j].insert(pair<int,int>(force_data.face_type[k],0));
   }

   // Forces are computed only on boundary faces
   for(unsigned int i=0; i<grid.n_boundary_face; ++i)
   {
      int face_type = grid.face[i].type;

      for(unsigned int j=0; j<force.size(); ++j)
      {
         ForceData& force_data = param.force_data[j];
         for(unsigned int k=0; k<force_data.face_type.size(); ++k)
            if(force_data.face_type[k] == face_type)
            {
               force[j].face.push_back (i);
               ++nface[j][face_type];
            }
      }
   }

   // Check for mistakes
   bool ok = true;

   // Check that all forces have faces
   for(unsigned int i=0; i<force.size(); ++i)
      if(force[i].face.size() == 0)
      {
         cout << "Force " << param.force_data[i].name << " does not have any faces\n";
         ok = false;
      }

   // Check that face types actually were found in the grid
   for(unsigned int j=0; j<force.size(); ++j)
   {
      ForceData& force_data = param.force_data[j];
      for(unsigned int k=0; k<force_data.face_type.size(); ++k)
         if(nface[j][force_data.face_type[k]] == 0)
         {
            cout << "Force: " << param.force_data[j].name << " has face type ";
            cout << force_data.face_type[k] << " but it is not present in the grid\n";
            ok = false;
         }
   }

   if(!ok)
      abort ();
}

//------------------------------------------------------------------------------
// Compute forces
// TODO: Currently computes only pressure forces
//------------------------------------------------------------------------------
void FiniteVolume::compute_forces (const unsigned int iter)
{
   force_file << setw(6) << iter << " " << scientific << setw(15);

   for(unsigned int i=0; i<force.size(); ++i)
   {
      force[i].value = 0.0;

      for(unsigned int j=0; j<force[i].face.size(); ++j)
      {
         const unsigned int face_no = force[i].face[j];
         vector<PrimVar> state(2);
         reconstruct ( face_no, false, state );
         force[i].value += grid.face[face_no].normal * state[0].pressure;
      }

      force_file << force[i].value.x << " " 
                 << force[i].value.y << " " 
                 << force[i].value.z << " ";
   }

   force_file << endl;
}
