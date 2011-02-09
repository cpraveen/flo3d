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
// Read parameters from file TODO
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
   assert (time_scheme == "rk1" || time_scheme == "rk3");
   if(time_scheme=="rk1") n_rks = 1;
   if(time_scheme=="rk3") n_rks = 3;

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
   checkString (input, "}");
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
   assert (material.gamma > 0.0);

   skipComment (fin);
   fin >> input;
   checkString (input, "gas_const");
   fin >> material.gas_const;
   assert (material.gas_const > 0.0);

   skipComment (fin);
   fin >> input;
   checkString (input, "model");
   fin >> material.model;
   assert (material.model == "euler");

   skipComment (fin);
   fin >> input;
   checkString (input, "flux");
   fin >> material.flux;
   assert (material.flux == "roe");

   skipComment (fin);
   fin >> input;
   checkString (input, "}");
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
   fin >> prim_inf.density;
   assert (prim_inf.density > 0.0);

   skipComment (fin);
   fin >> input;
   checkString (input, "velocity");
   fin >> prim_inf.velocity.x;
   fin >> prim_inf.velocity.y;
   fin >> prim_inf.velocity.z;

   skipComment (fin);
   fin >> input;
   checkString (input, "pressure");
   fin >> prim_inf.pressure;
   assert (prim_inf.pressure > 0.0);

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

   while (input!="}")
   {
      int b_type;
      fin >> b_type >> input;
      assert (b_type != -1); // -1 used for interior faces

      skipComment (fin);

      if(input=="farfield")
         bc_type.insert (pair<int,BCType>(b_type, farfield));
      else if(input=="slip")
         bc_type.insert (pair<int,BCType>(b_type, slip));
      else
      {
         cout << "   Unknown boundary type " << input << endl;
         abort ();
      }

      // Check if end of section is reached
      streampos curr_pos = fin.tellg ();
      fin >> input;
      if(input != "}")
         fin.seekg(curr_pos);
   }

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

   while (input != "}")
   {
      fin >> input;
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
   checkString (input, "variables");

   fin >> input;
   checkString (input, "{");

   fin >> input;
   while (input!="}")
   {
      assert (input=="density" || input=="velocity" || input=="pressure" ||
              input=="mach");
      write_variables.push_back (input);

      // Check if end of section is reached
      streampos curr_pos = fin.tellg ();
      fin >> input;
      if(input != "}")
      {
         fin.seekg(curr_pos);
         fin >> input;
      }
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
