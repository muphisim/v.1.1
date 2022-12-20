//+
lsca = 0.4;
c = 1;
L = 3;
Point(1) = {0, -c, 0, lsca};
Point(2) = {0.25*L, -c, 0, lsca};
Point(3) = {0.75*L, -c, 0, lsca};
Point(4) = {L, -c, 0, lsca};
Point(5) = {0, c, 0, lsca};
Point(6) = {0.25*L, c, 0, lsca};
Point(7) = {0.75*L, c, 0, lsca};
Point(8) = {L, c, 0, lsca};

//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 8};
//+
Line(5) = {8, 7};
//+
Line(6) = {7, 6};
//+
Line(7) = {6, 5};
//+
Line(8) = {5, 1};
//+
Line(9) = {6, 2};
//+
Line(10) = {7, 3};
//+
Curve Loop(1) = {7, 8, 1, -9};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, 2, -10, 6};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {10, 3, 4, 5};
//+
Plane Surface(3) = {3};
//+
Physical Surface(11) = {2};
//+
Physical Surface(10) = {3, 1};
//+
Physical Curve(12) = {8};
//+
Physical Curve(13) = {4};
