L = 1.0;  // plate length
s = 0.01; // semi span
H = 0.5; // height

n1 = 50;  // normal to BL
n2 = 50;  // along plate
p  = 1.1; // geometric progression

m  = 6;   // no. of layers in third direction

// First cell height
dh = H*(p - 1)/(p^n1 - 1);
Printf("First cell height = %e\n", dh);

Point(1) = {0,0,-s};
Point(2) = {L,0,-s};
Point(3) = {L,H,-s};
Point(4) = {0,H,-s};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Ruled Surface(1) = {1};
Transfinite Surface(1) = {1,2,3,4};

Transfinite Line{2,-4} = n1 Using Progression p;
Transfinite Line{1,-3} = n2;

out[] = Extrude{0,0,2*s}
{
  Surface{1};
  Layers{m};
};

Physical Volume(100000)  = {1};

Physical Surface(100001) = {25};   // inlet
Physical Surface(100002) = {13};   // plate
Physical Surface(100003) = {21};   // top
Physical Surface(100004) = {1,26}; // sides
Physical Surface(100005) = {17};   // outlet

