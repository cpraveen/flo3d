
// length (x)
L = 0.4;

// height (z)
H = 0.05;

// coefficient of increase in height (z direction)
Coef = 1.10 ;

// Thickness (y)
t = 0.05 ;

// plate length (x)
PL = 0.3;

// density of nodes
d_in = 0.2;
d_out = 0.2;
d_plate = 0.2;

// number of nodes
n1 = 10; 
n2 = 50; 
n3 = 30; 


//////////////////////////////////
// compute first step length

sum = 1;
For i In {1:n3}
sum = sum + Coef^i;
EndFor

H1 = H / sum;
Printf ("First layer thickness =%f",H1);

///////////////////////////////////////////
// first line

inlet = newp;
Point(inlet) = {0, 0, 0, d_in};

outlet = newp;
Point(outlet) = {L, 0, 0, d_out};

plate = newp;
Point(plate) = {L-PL, 0, 0, d_plate};

Line_inlet_plate = newreg;
Line(Line_inlet_plate) = {inlet,plate};

Line_plate_outlet = newreg;
Line(Line_plate_outlet) = {plate,outlet};

/////////////////////////////////////////////
// define nodes

// Put points on Line_inlet_plate
Transfinite Line{-Line_inlet_plate} = n1  Using Progression 1.15; 

// Put points on Line_plate_outlet
Transfinite Line{Line_plate_outlet} = n2  Using Progression 1.0; 

///////////////////////////////////////////
// extrude to mesh the plane

out[] = Extrude {0,t,0} { 
  Line{Line_inlet_plate,Line_plate_outlet}; Layers{6}; 
};

///////////////////////////////////////////
// extrude to mesh the volume

// first step

S1 = out[1];
S2 = out[5];
Hi = H1;

list_boundary_bottom[1] = S1;
list_boundary_plate[1] = S2;

Printf("Distance between nodes :");

For i In {1:n3}

tmp[] = Extrude {0,0,Hi} { 
  Surface{S1,S2}; Layers{1}; 
};

// compute the next step length

Printf ("%f",Hi);
Hi = Hi * Coef;

// attribute surfaces for the next step

S1 = tmp[0];
S2 = tmp[6];

// collect the physical entities

list_volume[2*i-1] = tmp[1];
list_volume[2*i] = tmp[7];

list_boundary_right[2*i-1] = tmp[2];
list_boundary_right[2*i] = tmp[8];

list_boundary_left[2*i-1] = tmp[4];
list_boundary_left[2*i] = tmp[10];

list_boundary_outlet[i] = tmp[9];
list_boundary_inlet[i] = tmp[5];

EndFor

list_boundary_top[1] = tmp[0];
list_boundary_top[2] = tmp[6];


///////////////////////////////////////////
// definition of the physical volume and boundaries


interior = 100000;
Physical Volume(interior) = list_volume[];

inlet = 100001;
Physical Surface(inlet) = list_boundary_inlet[];
outlet = 100002;
Physical Surface(outlet) = list_boundary_outlet[];

left = 100003;
Physical Surface(left) = list_boundary_left[];
right = 100004;
Physical Surface(right) = list_boundary_right[];

plate = 100005;
Physical Surface(plate) = list_boundary_plate[];

bottom = 100006;
Physical Surface(bottom) = list_boundary_bottom[];

top = 100007;
Physical Surface(top) = list_boundary_top[];
