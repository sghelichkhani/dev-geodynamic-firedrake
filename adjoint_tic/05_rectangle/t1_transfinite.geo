Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {1, 1, 0};
Point(4) = {0, 1, 0};

Line(5) = {1, 2};
Line(6) = {3, 2};
Line(7) = {3, 4};
Line(8) = {4, 1};


// Define a close loop
Curve Loop(1) = {8, 5, -6, 7};

// defining Plane surface here is essential 
Plane Surface(1) = {1};

// this is the part for Transfinite
Transfinite Curve{5, 6, 7, 8} = 100;
Transfinite Surface{1} = {1, 2, 3, 4};

Physical Curve(1) = {5};
Physical Curve(2) = {6};
Physical Curve(3) = {7};
Physical Curve(4) = {8};

//Mesh.Smoothing = 100;
Recombine Surface{1};

Physical Surface("domain") = {1};

