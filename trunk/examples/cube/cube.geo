
Lx=2.0;
Ly=1.0;
Lz=1.0;

h = 0.1;

Point(1) = {0, 0, 0, h};
Point(2) = {0, Ly, 0, h};
Point(3) = {0, Ly, Lz, h};
Point(4) = {0, 0, Lz, h};

Point(5) = {Lx, 0, 0, h};
Point(6) = {Lx, Ly, 0, h};
Point(7) = {Lx, Ly, Lz, h};
Point(8) = {Lx, 0, Lz, h};


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

Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};

Physical Volume(100000) = {1, 21};

Physical Surface(100001) = {1}; //Inflow plane
Physical Surface(100002) = {6}; //Outflow plane
Physical Surface(100003) = {2,3,4,5}; //side plane
