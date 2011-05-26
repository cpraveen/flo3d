// Pipe radius = 1
// Pipe length = 10 (-5, +5) along z axis

lc = 0.1;

Point(1) = {0.0, 0.0, -5.0, lc};
Point(2) = {1.0, 0.0, -5.0, lc};
Point(3) = {0.0, 1.0, -5.0, lc};
Point(4) = {-1.0, 0.0, -5.0, lc};
Point(5) = {0.0, -1.0, -5.0, lc};

Point(6) = {0.0, 0.0, 5.0, lc};
Point(7) = {1.0, 0.0, 5.0, lc};
Point(8) = {0.0, 1.0, 5.0, lc};
Point(9) = {-1.0, 0.0, 5.0, lc};
Point(10) = {0.0, -1.0, 5.0, lc};

Circle(11) = {2,1,3};
Circle(12) = {3,1,4};
Circle(13) = {4,1,5};
Circle(14) = {5,1,2};

Circle(15) = {8,6,7};
Circle(16) = {9,6,8};
Circle(17) = {10,6,9};
Circle(18) = {7,6,10};

Line(30) = {2, 7};
Line(31) = {3, 8};
Line(32) = {4, 9};
Line(33) = {5, 10};

Line(34) = {7, 2};
Line(35) = {8, 3};
Line(36) = {9, 4};
Line(37) = {10, 5};

Line Loop(38) = {34,11,-35,15};
Ruled Surface(39) = {38};
Line Loop(40) = {-36,16,35,12};
Ruled Surface(41) = {40};
Line Loop(42) = {36,13,-37,17};
Ruled Surface(43) = {42};
Line Loop(44) = {37,14,-34,18};
Ruled Surface(45) = {44};
Line Loop(46) = {16,15,18,17};
Plane Surface(47) = {46};
Line Loop(48) = {13,14,11,12};
Plane Surface(49) = {48};
Surface Loop(50) = {41,43,45,39,47,49};


Volume(51) = {50};

Physical Volume(52) = {51};
Physical Surface(100001) = {41,43,45,39}; // pipe wall
Physical Surface(100002) = {47};  // outlet
Physical Surface(100003) = {49};  // inlet

