#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "writer.h"

using namespace std;

// Add integer data
void Writer::attach_data (vector<unsigned int>& data, string name)
{
   integer_data.push_back (&data);
   integer_data_name.push_back (name);
}

// Write data to vtk file
void Writer::output_vtk (string filename)
{
   ofstream vtk;
   vtk.open (filename.c_str());
   vtk.close ();
}
