#include <iostream>
#include <cmath>
#include <map>
#include <cassert>
#include <cstdlib>
#include "parameter.h"

using namespace std;

//------------------------------------------------------------------------------
// Read parameters from file
//------------------------------------------------------------------------------
void Parameter::read ()
{
   cout << "Reading input file " << file << endl;
   Reader fin(file);

   read_grid (fin);
   read_numeric (fin);
   read_material (fin);
   read_initial_condition (fin);
   read_boundary (fin);
   read_integrals (fin);
   read_output (fin);
}

//------------------------------------------------------------------------------
// Read grid section
//------------------------------------------------------------------------------
void Parameter::read_grid (Reader &fin)
{
   cout << "  Reading grid section\n";

   string input;

   fin.begin_section ("grid");

   fin.entry ("type");
   fin >> input;
   if(input=="gmsh")
      grid_type = gmsh;
   else
   {
      cout << "   Unknown file type " << input << endl;
      abort ();
   }

   fin.entry ("file");
   fin >> grid_file;

   fin.end_section ();
}

//------------------------------------------------------------------------------
// Read numeric section
//------------------------------------------------------------------------------
void Parameter::read_numeric (Reader &fin)
{
   cout << "  Reading numeric section\n";

   string input;

   fin.begin_section ("numeric");

   fin.entry ("time_mode");
   fin >> time_mode;
   assert (time_mode == "steady" || time_mode == "unsteady");

   fin.entry ("time_scheme");
   fin >> time_scheme;
   assert (time_scheme == "rk1" || time_scheme == "rk3" || time_scheme == "lusgs");
   if(time_scheme=="rk1")   n_rks = 1;
   if(time_scheme=="rk3")   n_rks = 3;
   if(time_scheme=="lusgs") n_rks = 1; 

   fin.entry ("cfl");
   fin >> cfl;
   assert (cfl > 0.0);

   fin.entry ("max_iter");
   fin >> max_iter;
   assert (max_iter > 0);

   fin.entry ("final_time");
   fin >> final_time;
   assert (final_time > 0.0);

   fin.entry ("min_residue");
   fin >> min_residue;
   assert (min_residue > 0.0);

   fin.entry ("reconstruct");
   fin >> input;
   if(input == "first")
      reconstruct_scheme = Parameter::first;
   else if(input == "second")
      reconstruct_scheme = Parameter::second;
   else if(input == "secondF")
      reconstruct_scheme = Parameter::secondF;
   else if(input == "limited")
      reconstruct_scheme = Parameter::limited;
   else if(input == "limitedF")
      reconstruct_scheme = Parameter::limitedF;
   else if(input == "jameson")
   {
      reconstruct_scheme = Parameter::jameson;
      fin >> lim_power;
      assert (lim_power >= 1.0);
      assert (lim_power <= 3.0);
   }
   else
   {
      cout << "read_numeric: unknown reconstruction scheme " << input << endl;
      abort ();
   }

   fin.end_section ();

   // Some parameter checks
   if(time_scheme == "lusgs")
      assert (time_mode != "unsteady");
}

//------------------------------------------------------------------------------
// Read material section
//------------------------------------------------------------------------------
void Parameter::read_material (Reader &fin)
{
   cout << "  Reading material section\n";

   string input;

   fin.begin_section ("material");

   fin.entry ("gamma");
   fin >> material.gamma;
   assert (material.gamma > 1.0);

   fin.entry ("gas_const");
   fin >> material.gas_const;
   assert (material.gas_const > 0.0);

   fin.entry ("viscosity");
   fin >> input;
   if(input == "constant")
   {
      material.mu_model = Material::mu_constant;
      fin >> material.mu_ref;
      assert (material.mu_ref >= 0.0);
   }
   else if(input == "sutherland")
   {
      material.mu_model = Material::mu_sutherland;
      fin >> material.mu_ref;
      fin >> material.T_ref;
      fin >> material.T_0;
      assert (material.mu_ref > 0.0);
      assert (material.T_ref > 0.0);
      assert (material.T_0 > 0.0);
   }
   else
   {
      cout << "read_material: unknown viscosity type " << input << endl;
      abort ();
   }

   fin.entry ("prandtl");
   fin >> material.prandtl;
   assert (material.prandtl > 0.0);

   fin.entry ("model");
   fin >> input;
   if(input == "euler")
      material.model = Material::euler;
   else if(input == "ns")
      material.model = Material::ns;
   else
   {
      cout << "read_material: unknown flow model " << input << endl;
      abort ();
   }

   fin.entry ("flux");
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

   fin.end_section ();

   material.initialize ();
}

//------------------------------------------------------------------------------
// Read initial condition
//------------------------------------------------------------------------------
void Parameter::read_initial_condition (Reader &fin)
{
   cout << "  Reading initial condition section\n";

   string input;

   fin.begin_section ("initial_condition");

   fin.entry ("density");
   fin.getline (input);
   initial_condition.add ("density", input);

   fin.entry ("xvelocity");
   fin.getline (input);
   initial_condition.add ("xvelocity", input);

   fin.entry ("yvelocity");
   fin.getline (input);
   initial_condition.add ("yvelocity", input);

   fin.entry ("zvelocity");
   fin.getline (input);
   initial_condition.add ("zvelocity", input);

   fin.entry ("pressure");
   fin.getline (input);
   initial_condition.add ("pressure", input);

   fin.end_section ();
}

//------------------------------------------------------------------------------
// Read boundary conditions
//------------------------------------------------------------------------------
void Parameter::read_boundary (Reader &fin)
{
   cout << "  Reading boundary section\n";

   string input;

   fin.begin_section ("boundary");

   while (!fin.eos())
   {
      vector<int> f_type;
      while(!fin.bos())
      {
         int tmp;
         fin >> tmp;
         f_type.push_back (tmp);
      }

      fin.entry ("type");
      string bc_type;
      fin >> bc_type;

      vector<string> variable, function;
      while(!fin.eos())
      {
         fin >> input;
         variable.push_back (input);
         fin.getline(input);
         function.push_back (input);
      }
      BoundaryCondition bc(material, bc_type, variable, function);
      for(unsigned int i=0; i<f_type.size(); ++i)
         boundary_condition.insert (pair<int,BoundaryCondition>(f_type[i], bc));
   }
}

//------------------------------------------------------------------------------
// Read integrals section
//------------------------------------------------------------------------------
void Parameter::read_integrals (Reader &fin)
{
   cout << "  Reading integrals section\n";

   string input;

   fin.begin_section ("integrals");

   while (!fin.eos())
   {
      fin >> input;
      if(input=="force") // integral type is force
      {
         fin.entry ("{");

         force_data.resize( force_data.size() + 1 );
         fin >> force_data.back().name; // name to identify force

         while (!fin.eos())
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
void Parameter::read_output (Reader &fin)
{
   cout << "  Reading output section\n";

   string input;

   fin.begin_section ("output");

   fin.entry ("format");
   fin >> write_format;
   assert (write_format == "vtk");

   fin.entry ("frequency");
   fin >> write_frequency;
   assert (write_frequency > 0);

   fin.entry ("vertex");
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
   

   fin.begin_section ("variables");

   while (!fin.eos())
   {
      fin >> input;
      assert (input=="density" || input=="velocity" || input=="pressure" ||
              input=="mach");
      write_variables.push_back (input);
   }

   fin.entry ("restart");
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

   fin.end_section ();
}
