#include <iostream>
#include <cassert>
#include <cstring>

extern bool debug;

using namespace std;

void process_command_line (int   argc,
                           char* argv[],
                           int& ifile)
{
   assert (argc >= 3);

   // By default, debug is off
   debug = false;

   int i = 1;

   while (i < argc)
   {
      if(strcmp(argv[i],"-d")==0)
      {
         debug = true;
      }
      else if(strcmp(argv[i],"-i")==0)
      {
         assert (i+1 < argc); // check that there is another argument
         ifile = i+1;
         ++i;
      }
      else
      {
         cout << "Unknown command line flag: " << argv[i] << endl;
         cout << "Valid flags are:\n";
         cout << "   -i filename   Specify input file name\n";
         cout << "   -d            Enable debug mode\n";
         abort ();
      }

      ++i;
   }
}
