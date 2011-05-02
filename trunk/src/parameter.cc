#include <iostream>
#include <cmath>
#include <map>
#include <cassert>
#include <cstdlib>
#include "parameter.h"

using namespace std;

// Function declaration
ifstream& skipComment (ifstream &strm);

//------------------------------------------------------------------------------
// Check that two strings match
//------------------------------------------------------------------------------
void checkString (const string& str1, const string& str2)
{
   if(str1 != str2)
   {
      cout << "   Expecting " << str2 << " but found " << str1 << endl;
      abort ();
   }
}

//------------------------------------------------------------------------------
// Return true if end of section string "}" is encountered
// otherwise, put undo read and return false
//------------------------------------------------------------------------------
bool eos (ifstream& f)
{
   streampos curr_pos = f.tellg ();

   string input;
   f >> input;

   if(input != "}")
   {
      f.seekg (curr_pos);
      return false;
   }
   else
      return true;
}

//------------------------------------------------------------------------------
// Read parameters from file
//------------------------------------------------------------------------------
void Parameter::read ()
{
   cout << "Reading input file " << file << endl;
   fin.open (file);
   assert (fin.is_open());

   read_grid ();
   read_numeric ();
   read_material ();
   read_initial_condition ();
   read_boundary ();
   read_integrals ();
   read_output ();

   fin.close ();
}

//------------------------------------------------------------------------------
// Read grid section
//------------------------------------------------------------------------------
void Parameter::read_grid ()
{
   cout << "  Reading grid section\n";

   string input;

   skipComment (fin);
   fin >> input;
   checkString (input, "grid");

   fin >> input;
   checkString (input, "{");

   skipComment (fin);
   fin >> input;
   checkString (input, "type");
   fin >> input;
   if(input=="gmsh")
      grid_type = gmsh;
   else
   {
      cout << "   Unknown file type " << input << endl;
      abort ();
   }

   skipComment (fin);
   fin >> input;
   checkString (input, "file");
   fin >> grid_file;

   skipComment (fin);
   fin >> input;
   checkString (input, "}");
}

//------------------------------------------------------------------------------
// Read numeric section
//------------------------------------------------------------------------------
void Parameter::read_numeric ()
{
   cout << "  Reading numeric section\n";

   string input;

   skipComment (fin);
   fin >> input;
   checkString (input, "numeric");

   fin >> input;
   checkString (input, "{");

   skipComment (fin);
   fin >> input;
   checkString (input, "time_mode");
   fin >> time_mode;
   assert (time_mode == "steady" || time_mode == "unsteady");

   skipComment (fin);
   fin >> input;
   checkString (input, "time_scheme");
   fin >> time_scheme;
   assert (time_scheme == "rk1" || time_scheme == "rk3" || time_scheme == "lusgs");
   if(time_scheme=="rk1")   n_rks = 1;
   if(time_scheme=="rk3")   n_rks = 3;
   if(time_scheme=="lusgs") n_rks = 1; 

   skipComment (fin);
   fin >> input;
   checkString (input, "cfl");
   fin >> cfl;
   assert (cfl > 0.0);

   skipComment (fin);
   fin >> input;
   checkString (input, "max_iter");
   fin >> max_iter;
   assert (max_iter > 0);

   skipComment (fin);
   fin >> input;
   checkString (input, "final_time");
   fin >> final_time;
   assert (final_time > 0.0);

   skipComment (fin);
   fin >> input;
   checkString (input, "min_residue");
   fin >> min_residue;
   assert (min_residue > 0.0);

   skipComment (fin);
   fin >> input;
   checkString (input, "reconstruct");
   fin >> input;
   if(input == "first")
      reconstruct_scheme = Parameter::first;
   else if(input == "second")
      reconstruct_scheme = Parameter::second;
   else if(input == "limited")
      reconstruct_scheme = Parameter::limited;
   else
   {
      cout << "read_numeric: unknown reconstruction scheme " << input << endl;
      abort ();
   }

   skipComment (fin);
   fin >> input;
   checkString (input, "}");

   // Some parameter checks
   if(time_scheme == "lusgs")
      assert (time_mode != "unsteady");
}

//------------------------------------------------------------------------------
// Read material section
//------------------------------------------------------------------------------
void Parameter::read_material ()
{
   cout << "  Reading material section\n";

   string input;

   skipComment (fin);
   fin >> input;
   checkString (input, "material");

   fin >> input;
   checkString (input, "{");

   skipComment (fin);
   fin >> input;
   checkString (input, "gamma");
   fin >> material.gamma;
   assert (material.gamma > 1.0);

   skipComment (fin);
   fin >> input;
   checkString (input, "gas_const");
   fin >> material.gas_const;
   assert (material.gas_const > 0.0);

   skipComment (fin);
   fin >> input;
   checkString (input, "prandtl");
   fin >> material.prandtl;
   assert (material.prandtl > 0.0);

   skipComment (fin);
   fin >> input;
   checkString (input, "model");
   fin >> material.model;
   assert (material.model == "euler" || material.model == "ns");

   skipComment (fin);
   fin >> input;
   checkString (input, "flux");
   fin >> input;
   if(input == "roe")
      material.flux_scheme = Material::roe;
   else if(input == "kfvs")
      material.flux_scheme = Material::kfvs;
   else
   {
      cout << "read_material:: unknown flux scheme: " << input << endl;
      abort ();
   }

   skipComment (fin);
   fin >> input;
   checkString (input, "}");

   material.initialize ();
}

//------------------------------------------------------------------------------
// Read initial condition
//------------------------------------------------------------------------------
void Parameter::read_initial_condition ()
{
   cout << "  Reading initial condition section\n";

   string input;

   skipComment (fin);
   fin >> input;
   checkString (input, "initial_condition");

   fin >> input;
   checkString (input, "{");

   skipComment (fin);
   fin >> input;
   checkString (input, "density");
   getline (fin, input);
   initial_condition.add ("density", input);

   skipComment (fin);
   fin >> input;
   checkString (input, "xvelocity");
   getline (fin, input);
   initial_condition.add ("xvelocity", input);

   skipComment (fin);
   fin >> input;
   checkString (input, "yvelocity");
   getline (fin, input);
   initial_condition.add ("yvelocity", input);

   skipComment (fin);
   fin >> input;
   checkString (input, "zvelocity");
   getline (fin, input);
   initial_condition.add ("zvelocity", input);

   skipComment (fin);
   fin >> input;
   checkString (input, "pressure");
   getline (fin, input);
   initial_condition.add ("pressure", input);

   skipComment (fin);
   fin >> input;
   checkString (input, "}");
}

//------------------------------------------------------------------------------
// Read boundary conditions
//------------------------------------------------------------------------------
void Parameter::read_boundary ()
{
   cout << "  Reading boundary section\n";

   string input;

   skipComment (fin);
   fin >> input;
   checkString (input, "boundary");

   fin >> input;
   checkString (input, "{");

   skipComment (fin);

   while (!eos(fin))
   {
      int b_type;
      fin >> b_type >> input;
      assert (b_type != -1); // -1 used for interior faces


      if(input=="farfield")
         bc_type.insert (pair<int,BCType>(b_type, farfield));
      else if(input=="slip")
         bc_type.insert (pair<int,BCType>(b_type, slip));
      else if(input=="noslip")
         bc_type.insert (pair<int,BCType>(b_type, noslip));
      else if(input=="inlet")
         bc_type.insert (pair<int,BCType>(b_type, inlet));
      else if(input=="outlet")
         bc_type.insert (pair<int,BCType>(b_type, outlet));
      else if(input=="pressure")
         bc_type.insert (pair<int,BCType>(b_type, pressure));
      else
      {
         cout << "   Unknown boundary type " << input << endl;
         abort ();
      }
      
      // Read primitive state for this boundary
      PrimVar state;
      fin >> state.density
          >> state.velocity.x
          >> state.velocity.y
          >> state.velocity.z
          >> state.pressure;

      // Density and pressure must be positive
      if(input == "inlet" || input == "farfield")
         assert (state.density  > 0.0);
      if(input == "inlet" || input == "farfield" || input == "pressure")
         assert (state.pressure > 0.0);
      
      bc_state.insert (pair<int,PrimVar>(b_type, state));

      skipComment (fin);

   }

   // We added -1 as interior type
   bc_type.insert (pair<int,BCType>(-1, interior));

}

//------------------------------------------------------------------------------
// Read integrals section
//------------------------------------------------------------------------------
void Parameter::read_integrals ()
{
   cout << "  Reading integrals section\n";

   string input;

   skipComment (fin);
   fin >> input;
   checkString (input, "integrals");

   fin >> input;
   checkString (input, "{");

   while (!eos(fin))
   {
      fin >> input;
      if(input=="force") // integral type is force
      {
         fin >> input;
         checkString (input, "{");

         force_data.resize( force_data.size() + 1 );
         fin >> force_data.back().name; // name to identify force

         while (!eos(fin))
         {
            int face_type;
            fin >> face_type;
            force_data.back().face_type.push_back(face_type);
         }

         // Check that there was atleast one face type given
         assert (force_data.back().face_type.size() > 0);
      }
   }

}

//------------------------------------------------------------------------------
// Read output section
//------------------------------------------------------------------------------
void Parameter::read_output ()
{
   cout << "  Reading output section\n";

   string input;

   skipComment (fin);
   fin >> input;
   checkString (input, "output");

   fin >> input;
   checkString (input, "{");

   skipComment (fin);
   fin >> input;
   checkString (input, "format");
   fin >> write_format;
   assert (write_format == "vtk");

   skipComment (fin);
   fin >> input;
   checkString (input, "frequency");
   fin >> write_frequency;
   assert (write_frequency > 0);

   skipComment (fin);
   fin >> input;
   checkString (input, "vertex");
   fin >> input;
   if(input == "true")
      write_vertex_variables = true;
   else if(input == "false")
      write_vertex_variables = false;
   else
   {
      cout << "read_output: unknown option = " << input << endl;
      abort ();
   }
   

   skipComment (fin);
   fin >> input;
   checkString (input, "variables");

   fin >> input;
   checkString (input, "{");

   while (!eos(fin))
   {
      fin >> input;
      assert (input=="density" || input=="velocity" || input=="pressure" ||
              input=="mach");
      write_variables.push_back (input);
   }

   skipComment (fin);
   fin >> input;
   checkString (input, "restart");
   fin >> input;
   if(input=="false")
      write_restart = false;
   else if(input=="true")
      write_restart = true;
   else
   {
      cout << "   Unknown input: " << input << endl;
      abort ();
   }

   skipComment (fin);
   fin >> input;
   checkString (input, "}");
}
