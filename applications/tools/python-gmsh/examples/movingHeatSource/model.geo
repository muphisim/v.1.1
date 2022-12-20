unit = 1.e-3;
lsca = 0.025*unit;
eleSize = 0.005*unit;
Point(1) = {0,0,0,lsca};
Point(2) = {0.2*unit,0.,0,lsca};
//+
Line(1) = {1, 2};
//+
Extrude {0, 0.3*unit, 0} {
  Curve{1}; Layers {Ceil(0.3*unit/eleSize)}; Recombine;
}
//+
Extrude {0, 1*unit, 0} {
  Curve{2}; Layers {Ceil(1*unit/eleSize)}; Recombine;
}
//+
Extrude {0, 0.3*unit, 0} {
  Curve{6}; Layers {Ceil(0.3*unit/eleSize)}; Recombine;
}
//+
Extrude {0, 0, 0.2*unit} {
  Surface{13}; Surface{9}; Surface{5}; Layers {Ceil(0.2*unit/eleSize)}; Recombine;
}
//+
Physical Volume(80) = {3};
//+
Physical Volume(81) = {2};
//+
Physical Volume(82) = {1};
//+
Physical Surface(83) = {66};
//+
Physical Surface(84) = {30};
//+
Transfinite Curve {1} =15 Using Progression 1;
//+
Physical Surface(85) = {30, 35, 34, 26, 48, 57, 56, 78, 79, 70, 66};
//+
Physical Surface(86) = {5, 9, 13};
