r1 = 0.1; // radius of inner cylinder
r2 = 0.2; // radius of outer cylinder
h  = 0.1; // height of cylinder

cl = 0.01;

Point(1) = {r1, 0, 0, cl};
Point(2) = {r2, 0, 0, cl};
Point(3) = {r2, h, 0, cl};
Point(4) = {r1, h, 0, cl};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Extrude{ {0,1,0}, {0,0,0}, Pi/2}
{
   Surface{1};
}

Symmetry {0,0,1,0}{ Duplicata{Volume{1};}}
Symmetry {1,0,0,0}{ Duplicata{Volume{1};}}
Symmetry {1,0,0,0}{ Duplicata{Volume{27};}}

Physical Volume(1000) = {1, 27, 54, 81};
Physical Surface(1001) = {13, 65, 92, 38}; //bottom
Physical Surface(1002) = {21,75,102, 48}; //top
Physical Surface(1003) = {25, 80, 53, 107}; //inner cylinder
Physical Surface(1004) = {43, 97, 17, 70}; //outer cylinder
