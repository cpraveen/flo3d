r1 = 0.1; // radius of inner cylinder
r2 = 0.2; // radius of outer cylinder
h  = 0.1; // height of cylinder

n1     = 25; // points along radial
n2     = 25; // points along axial
ntheta = 25; // layers along tangential in one quadrant

p      = 0.1;// progression

Point(1) = {r1, 0, 0};
Point(2) = {r2, 0, 0};
Point(3) = {r2, h, 0};
Point(4) = {r1, h, 0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Ruled Surface(1) = {1};
Transfinite Surface(1) = {1,2,3,4};
Transfinite Line{1,3} = n1 Using Bump p;
Transfinite Line{2,4} = n2 Using Bump p;

Extrude{ {0,1,0}, {0,0,0}, Pi/2}
{
   Surface{1}; Layers{ntheta};
}

Extrude{ {0,1,0}, {0,0,0}, Pi/2}
{
   Surface{26}; Layers{ntheta};
}

Extrude{ {0,1,0}, {0,0,0}, Pi/2}
{
   Surface{48}; Layers{ntheta};
}

Extrude{ {0,1,0}, {0,0,0}, Pi/2}
{
   Surface{70}; Layers{ntheta};
}

Physical Volume (1000) = {1,2,3,4};

Physical Surface(1001) = {57,35,79,13}; // bottom
Physical Surface(1002) = {65,43,87,21}; // top
Physical Surface(1003) = {25,47,69,91}; // inner cylinder
Physical Surface(1004) = {39,17,83,61}; // outer cylinder
