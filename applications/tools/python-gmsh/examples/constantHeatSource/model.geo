unit = 1.;
lsca = 0.1;
eleSize = 0.01;
Point(1) = {0,0,0,lsca};
Point(2) = {0.1,0.,0,lsca};
//+
Line(1) = {1, 2};
//+
Extrude {0, 0.6, 0} {
  Curve{1}; Layers {Ceil(0.6/eleSize)}; Recombine;
}
//+
Extrude {0, 0.8, 0} {
  Curve{2}; Layers {Ceil(0.8/eleSize)}; Recombine;
}
//+
Extrude {0, 0.6, 0} {
  Curve{6}; Layers {Ceil(0.6/eleSize)}; Recombine;
}
//+
Extrude {0, 0, 0.1} {
  Surface{13}; Surface{9}; Surface{5}; Layers {4}; Recombine;
}
//+
Physical Volume(80) = {1,2,3};
//+
Physical Surface(83) = {66};
//+
Physical Surface(84) = {30};
//+
Transfinite Curve {1} =Ceil(1/eleSize) Using Progression 1;
