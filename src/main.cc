#include <iostream>
#include <cstdlib>
#include <cassert>
#include "fv.h"

using namespace std;

int main(int argc, char* argv[])
{
   assert (argc==2);
   FiniteVolume problem (argv[1]);
   problem.run ();
}
