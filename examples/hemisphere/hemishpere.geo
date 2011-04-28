radius = 2.1;
cellSize=0.05;
cellSize1 = 0.05;  // inner hemisphere
pio2=Pi/2;
radius1 = 1.0;

// create inner 1/8 shell
Point(41) = {0, 0, 0, cellSize};
Point(42) = {-radius, 0, 0, cellSize};
Point(43) = {0, radius, 0, cellSize};
Point(46) = {0, 0, radius, cellSize};
Circle(21) = {42, 41, 43};
Circle(25)= {46,41,42};
Circle(26)= {46,41,43};
Line Loop(10) = {21, -26, 25} ;
Ruled Surface (60) = {10};


Point(411) = {0, 0, 0, cellSize1};
Point(421) = {-radius1, 0, 0, cellSize1};
Point(431) = {0, radius1, 0, cellSize1};
Point(461) = {0, 0, radius1, cellSize1};
Circle(211) = {421, 411, 431};
Circle(251)= {461,411,421};
Circle(261)= {461,411,431};
Line Loop(101) = {211, -261, 251} ;
Ruled Surface (601) = {101};


t1[] = Rotate {{0,0,1},{0,0,0},pio2} {Duplicata{Surface{60};}};
t2[] = Rotate {{0,0,1},{0,0,0},pio2*2} {Duplicata{Surface{60};}};
t3[] = Rotate {{0,0,1},{0,0,0},pio2*3} {Duplicata{Surface{60};}};


t11[] = Rotate {{0,0,1},{0,0,0},pio2} {Duplicata{Surface{601};}};
t21[] = Rotate {{0,0,1},{0,0,0},pio2*2} {Duplicata{Surface{601};}};
t31[] = Rotate {{0,0,1},{0,0,0},pio2*3} {Duplicata{Surface{601};}};





// create entire inner and outer shell
Surface Loop(100)={60,t1[0],t2[0],t3[0]};
Surface Loop(200)={601,t11[0],t21[0],t31[0]};


Line(1000) = {42,421};
Line(1001) = {43,431};

Line Loop(10001) = {1000,211 , -1001, -21 };

Plane Surface (60001) = {10001};

Line(1002) = {462,464};

Line(1003) = {465,463};

Line Loop(10002) = {1000,-613,-1002,603};
Plane Surface (60002) = {10002};


Line Loop(10003) = {1002,-617,1003,607};
Plane Surface (60003) = {10003};


Line Loop(10004) = {1003,-611,1001,621};
Plane Surface (60004) = {10004};


Surface Loop(10000)={60,601,t1[0],t2[0],t3[0],t11[0],t21[0],t31[0],60001,60002,60003,60004};



//create volume between shells
Volume(100000)={10000};
// 60001,60002,60003,60004 are outlet surface
// 60,t1,t2,t3 are inlet
// 601, t11,t21,t31 is desired hemisphere 

Physical Surface(600001) = {60001,60002,60003,60004};  // outlet

Physical Surface(600002) = {60,t1[0],t2[0],t3[0] };
Physical Surface(600003) = {601,t11[0],t21[0],t31[0] };
Physical Volume( 10000001) = { 100000};



