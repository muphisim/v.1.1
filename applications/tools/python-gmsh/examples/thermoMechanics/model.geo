//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0.1, 0, 0, 1.0};
//+
Point(3) = {1, 0, 0, 1.0};
//+
Point(4) = {1, 1, 0, 1.0};
//+
Point(5) = {0, 1, 0, 1.0};
//+
Point(6) = {0, 0.1, 0, 1.0};
//+
Line(1) = {6, 5};
//+
Line(2) = {5, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 2};
//+
Circle(5) = {2, 1, 6};
//+
Curve Loop(1) = {1, 2, 3, 4, 5};
//+
Plane Surface(1) = {1};
//+
Physical Surface(6) = {1};
//+
Physical Curve(7) = {4};
//+
Physical Curve(8) = {3};
//+
Physical Curve(9) = {2};
//+
Physical Curve(10) = {1};
//+
Physical Curve(11) = {5};
//+
//+
MeshSize {6, 2} = 0.001;
//+
MeshSize {5, 4, 3} = 0.05;
