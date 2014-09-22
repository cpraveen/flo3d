// inner sphere radius is 1
R = 20;        //outer sphere radius
lc1 = 2*Pi/200; // characteristic length
lc2 = 2*Pi*R/20; // characteristic length

Point(1) = {0.0,0.0,0.0};
Point(2) = {1,0.0,0.0};
Point(3) = {0,1,0.0};
Circle(1) = {2,1,3};
Point(4) = {-1,0,0.0};
Point(5) = {0,-1,0.0};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Point(6) = {0,0,-1};
Point(7) = {0,0,1};
Circle(5) = {3,1,6};
Circle(6) = {6,1,5};
Circle(7) = {5,1,7};
Circle(8) = {7,1,3};
Circle(9) = {2,1,7};
Circle(10) = {7,1,4};
Circle(11) = {4,1,6};
Circle(12) = {6,1,2};
Line Loop(13) = {2,8,-10};
Ruled Surface(14) = {13};
Line Loop(15) = {10,3,7};
Ruled Surface(16) = {15};
Line Loop(17) = {-8,-9,1};
Ruled Surface(18) = {17};
Line Loop(19) = {-11,-2,5};
Ruled Surface(20) = {19};
Line Loop(21) = {-5,-12,-1};
Ruled Surface(22) = {21};
Line Loop(23) = {-3,11,6};
Ruled Surface(24) = {23};
Line Loop(25) = {-7,4,9};
Ruled Surface(26) = {25};
Line Loop(27) = {-4,12,-6};
Ruled Surface(28) = {27};
Surface Loop(29) = {28,26,16,14,20,24,22,18};

// Scale sphere surface to obtain outer boundary
out[] = Dilate {{0,0,0},R}
{
   Duplicata{ Surface{14,16,18,20,22,24,26,28}; }
};

// Outer boundary
Surface Loop(300) = {30,34,38,54,42,50,46,58};

Volume(31) = {300,29};

// Attractor to sphere surface
Field[1] = Attractor;
Field[1].FacesList = {14,16,18,20,22,24,26,28};

Field[2] = MathEval;
Field[2].F = Sprintf("(1-(F1/%g)^1.5)*%g + (F1/%g)^1.5*%g", R, lc1, R, lc2);
//Field[2].F = Sprintf("(1+F1^2)*%g", lc);

Background Field = 2;

// try also netgen:
// Mesh.Algorithm3D = 4;

Physical Volume(100000) = {31};
Physical Surface(100001) = {14,16,18,20,22,24,26,28}; // sphere surface
Physical Surface(100002) = {30,34,38,54,42,50,46,58}; // outer boundary
