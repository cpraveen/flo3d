L = 4.0; // extrusion in third direction
R = 2.0; // radius of outer cylinder
h = 0.1; // point density at inner cylinder
l1= 1*R;
l2=1*R;
C= Sqrt(R^2-0.5^2);
H=1*R;

Point(1) = {0,0,0,h};

Point(2) = {1,0,0,h};
Point(3) = {1+l2,0,0,h};
Point(4) = {1+l2,H,0,h};
Point(5) = {-l1,H,0,h};

Point(6) = {-l1,0,0,h};
Point(7) = {0.5,-C,h};

Circle(1) = {1,7,2};

Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};


Line Loop(1) = {1,2,3,4,5,6};

Plane Surface(2) = {1};



out[]=Extrude {0,0,.5} {
  Surface{2};
};

Physical Volume(100000) = {out[1]};  // volume
Physical Surface(100001) = {out[0]}; //outflow surface
Physical Surface(100002) = {out[2]}; //lower surface

Physical Surface(100003) = {out[3]}; //lower surface
Physical Surface(100004) = {out[4]}; //side surface

Physical Surface(100005) = {out[5]}; //upper surface
Physical Surface(100006) = {out[6]}; //side surface

Physical Surface(100007) = {out[7]}; //lower surface
Physical Surface(100008) = {2}; // inflow face

//printf( "out[0]");
//show { L};

  Printf("  %g   %g  %g    %g    %g  %g    %g    %g  " ,
	 out[0],out[1],out[2],out[3],out[4],out[5],out[6],out[7]) ;


/*Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};

Physical Volume(100000) = {1, 21};

Physical Surface(100001) = {1}; //Inflow plane
Physical Surface(100002) = {6}; //Outflow plane
Physical Surface(100003) = {2,3,4,5}; //side plane
*/




