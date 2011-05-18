// To run this:
//   gmsh -3 -optimize_netgen shocktube.geo
//

Lx=1;
Ly=.1;
Lz=.1;
Lx1=0.5; // Location of initial discontinuity

h = 0.01;

Point(1) = {0, 0, 0, h};
Point(2) = {0, Ly, 0, h};
Point(3) = {0, Ly, Lz, h};
Point(4) = {0, 0, Lz, h};

Point(5) = {Lx, 0, 0, h};
Point(6) = {Lx, Ly, 0, h};
Point(7) = {Lx, Ly, Lz, h};
Point(8) = {Lx, 0, Lz, h};

Point(9) = {Lx1, 0, 0, h};
Point(10) = {Lx1, Ly, 0, h};
Point(11) = {Lx1, Ly, Lz, h};
Point(12) = {Lx1, 0, Lz, h};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};


Line(5) = {1,5};
Line(6) = {5,6};
Line(7) = {6,2};

Line(8) = {6,7};
Line(9) = {7,3};

Line(10) = {4,8};
Line(11) = {8,7};

Line(12) = {5,8};

Line(13) = {9,10};
Line(14) = {10,11};
Line(15) = {11,12};
Line(16) = {12,9};


Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Line Loop(2) = {5,6,7,-1};
Plane Surface(2) = {2};

Line Loop(3) = {-7,8,9,-2};
Plane Surface(3) = {3};

Line Loop(4) = {10,11,9,3};
Plane Surface(4) = {4};

Line Loop(5) = {5,12,-10,4};
Plane Surface(5) = {5};

Line Loop(6) = {6,8,-11,-12};
Plane Surface(6) = {6};

Line Loop(7) = {13,14,15,16};
Plane Surface(7) = {7};

Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};

Symmetry{0,1,0,0}
{
  Duplicata{Surface{7}; Volume{1};}
}

Symmetry{0,0,1,0}
{
  Duplicata{Surface{7}; Volume{1};}
}

Symmetry{0,1,0,0}
{
  Duplicata{Surface{49}; Volume{54};}
}

Physical Volume(100000) = {1, 22, 54, 86};

Physical Surface(100001) = {1,23,55,87}; //Inflow plane
Physical Surface(100002) = {6,48,80,112}; //Outflow plane
Physical Surface(100003) = {4,38,3,65,33,97,70,102}; //side plane
