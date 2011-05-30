#include <iostream>
#include <cassert>
#include <cstring>
#include <cstdlib>

extern bool debug;
extern bool restart;
extern bool preprocess;
extern bool bounds;

using namespace std;

void show_options ();

//------------------------------------------------------------------------------
// Get command line flags and input file
//------------------------------------------------------------------------------
void process_command_line (int   argc,
                           char* argv[],
                           int&  ifile)
{ 
   if(argc < 3)
      show_options ();

   // By default, all are off
   debug      = false;
   restart    = false;
   preprocess = false;
   bounds     = false;

   int i = 1;
   bool found_input_file = false;

   while (i < argc)
   {
      if(strcmp(argv[i],"-d")==0)
      {
         debug = true;
      }
      else if(strcmp(argv[i],"-r")==0)
      {
         restart = true;
      }
      else if(strcmp(argv[i],"-p")==0)
      {
         preprocess = true;
      }
      else if(strcmp(argv[i],"-b")==0)
      {
         bounds = true;
      }
      else if(strcmp(argv[i],"-i")==0)
      {
         assert (i+1 < argc); // check that there is another argument
         ifile = i+1;
         ++i;
         found_input_file = true;
      }
      else
      {
         cout << "Unknown command line flag: " << argv[i] << endl;
         show_options ();
      }

      ++i;
   }

   if(!found_input_file)
      show_options ();
}

//------------------------------------------------------------------------------
// Print command line options available
//------------------------------------------------------------------------------
void show_options ()
{
   cout << "Valid flags are:\n";
   cout << "   -i filename   Specify input file name (required)\n";
   cout << "   -d            Enable debug mode (optional)\n";
   cout << "   -r            Read restart file for initial condition (optional)\n";
   cout << "   -p            Do everything but do not solve (optional)\n";
   cout << "   -b            Compute min/max range of solution (optional)\n";
   abort ();
}
