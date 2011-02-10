#include <iostream>
#include <cstdlib>
#include <cassert>
#include "fv.h"

bool debug;
bool restart;

using namespace std;

void process_command_line (int argc, char* argv[], int& ifile);

int main(int argc, char* argv[])
{
   cout << "Starting flo3d\n";   
   int ifile;
   process_command_line (argc, argv, ifile);

   FiniteVolume problem (argv[ifile]);
   problem.run ();
}
