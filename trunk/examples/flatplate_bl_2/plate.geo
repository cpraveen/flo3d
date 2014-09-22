l = 0.5;  // initial length
L = 1.0;  // plate length
s = 0.01; // semi span
H = 0.5; // height

n1 = 50;  // normal to BL
n2 = 50;  // along plate
n3 = 25;  // initial portion
p  = 1.1; // geometric progression

m  = 6;   // no. of layers in third direction

// First cell height
dh = H*(p - 1)/(p^n1 - 1);
Printf("First cell height = %e\n", dh);

Point(1) = {0,0,-s};
Point(2) = {L,0,-s};
Point(3) = {L,H,-s};
Point(4) = {0,H,-s};

Point(5) = {-l,0,-s};
Point(6) = {-l,H,-s};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {1,5};
Line(6) = {5,6};
Line(7) = {6,4};

Line Loop(1) = {1,2,3,4};
Ruled Surface(1) = {1};
Transfinite Surface(1) = {1,2,3,4};

Line Loop(2) = {-4,-7,-6,-5};
Ruled Surface(2) = {2};
Transfinite Surface(2) = {1,5,6,4};

Transfinite Line{2,-4,6} = n1 Using Progression p;
Transfinite Line{1,-3} = n2 Using Progression 1.05;
Transfinite Line{5,-7} = n3 Using Progression 1.1;

out[] = Extrude{0,0,2*s}
{
  Surface{1,2};
  Layers{m};
};

Physical Volume(100000)  = {1,2};

Physical Surface(100001) = {46};        // inlet
Physical Surface(100002) = {16};        // plate
Physical Surface(100003) = {42,24};     // top
Physical Surface(100004) = {1,2,29,51}; // sides
Physical Surface(100005) = {20};        // outlet
Physical Surface(100006) = {50};        // inlet mirror face on bottom

