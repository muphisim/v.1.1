unit = 1.;
lsca = 0.1;
eleSize = 0.025;
Point(1) = {0,0,0,lsca};
Point(2) = {0.05,0.,0,lsca};
//+
Line(1) = {1, 2};
//+
Extrude {0, 0.3, 0} {
  Curve{1}; Layers {Ceil(0.3/eleSize)}; Recombine;
}
//+
Extrude {0, 0.1, 0} {
  Curve{2}; Layers {Ceil(0.1/eleSize)}; Recombine;
}
//+
Extrude {0, 0.3, 0} {
  Curve{6}; Layers {Ceil(0.3/eleSize)}; Recombine;
}
//+
Extrude {0, 0, 0.05} {
  Surface{13}; Surface{9}; Surface{5}; Layers {4}; Recombine;
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
Transfinite Curve {1} =5 Using Progression 1;
