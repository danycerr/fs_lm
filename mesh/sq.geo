//+
l=1;
ndiv=8;
Point(1) = {0, 0, -0, 1.0};
Point(2) = {0, l, -0, 1.0};
Point(3) = {l, l, -0, 1.0};
Point(4) = {l, 0, -0, 1.0};
Point(5) = {l/2, l, -0, 1.0};
Point(6) = {l/2, 0, -0, 1.0};
//+
Line(1) = {1, 6};
//+
Line(2) = {6, 5};
//+
Line(3) = {5, 2};
//+
Line(4) = {2, 1};
//+
Line(5) = {6, 4};
//+
Line(6) = {5, 3};
//+
Line(7) = {3, 4};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {6, 7, -5, 2};

Plane Surface(2) = {2};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Curve {4, 2, 7} = ndiv +1 Using Progression 1;
//+
Transfinite Curve {3, 6, 5, 1} = ndiv/2 +1 Using Progression 1;
//+
Transfinite Surface {1} Alternated;
//+
Recombine Surface {1, 2};
//+
Physical Surface("1") = {1};
//+
Physical Surface("2") = {2};
//+
Physical Curve("10") = {2};
