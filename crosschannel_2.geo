Mesh.MshFileVersion = 2.2;
SetFactory("OpenCASCADE");

// ======================================================
// 1. CROSS GEOMETRY
// ======================================================

L = channel_length;
W = channel_height;

// Outer boundary points
Point(1) = {-L, -W/2, 0, elemsize};
Point(2) = {-W/2, -W/2, 0, elemsize};
Point(3) = {-W/2, -L, 0, elemsize};
Point(4) = { W/2, -L, 0, elemsize};
Point(5) = { W/2, -W/2, 0, elemsize};
Point(6) = { L, -W/2, 0, elemsize};
Point(7) = { L, W/2, 0, elemsize};
Point(8) = { W/2, W/2, 0, elemsize};
Point(9) = { W/2, L, 0, elemsize};
Point(10) = {-W/2, L, 0, elemsize};
Point(11) = {-W/2, W/2, 0, elemsize};
Point(12) = {-L, W/2, 0, elemsize};
Point(13) = {0, 0.0001, 0, 1.0};
Point(14) = {0, -0.0001, 0, 1.0};

// Lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,1};

Line Loop(100) = {1,2,3,4,5,6,7,8,9,10,11,12};
Plane Surface(200) = {100};

// ======================================================
// 2. CENTRAL SOLID DISK
// ======================================================
particle_y_offset = 0;
particle_x_offset = -3*1e-4;
Disk(300) = {particle_x_offset,particle_y_offset,0,cylinder_radius};

// ======================================================
// 3. FLUID DOMAIN
// ======================================================
BooleanDifference(400) = { Surface{200}; Delete; }{ Surface{300}; };
A = newp; Point(A) = { particle_x_offset, particle_y_offset, 0 };
B = newp; Point(B) = { -cylinder_radius + particle_x_offset, particle_y_offset, 0 };

_() = BooleanFragments{ Surface{300,400}; Point{A,B}; Delete; }{};

// ======================================================
// 4. BOUNDARY CLASSIFICATION
// ======================================================

bnd_fluid() = Abs(Boundary{ Surface{400}; });
bnd_solid() = Abs(Boundary{ Surface{300}; });

bnd_fluid -= bnd_solid();

inlet1() = {};
inlet2() = {};
outlet1() = {};
outlet2() = {};
wall() = {};

For i In {0:#bnd_fluid()-1}

pts[] = Boundary{ Curve{ bnd_fluid(i) }; };

xc = 0;
yc = 0;

For j In {0:#pts()-1}
coord[] = Point{ pts(j) };
xc += coord[0];
yc += coord[1];
EndFor

xc /= #pts();
yc /= #pts();

tol = 1e-12;

If (Abs(xc + L) < tol)
inlet1() += { bnd_fluid(i) };
ElseIf (Abs(xc - L) < tol)
inlet2() += { bnd_fluid(i) };
ElseIf (Abs(yc + L) < tol)
outlet1() += { bnd_fluid(i) };
ElseIf (Abs(yc - L) < tol)
outlet2() += { bnd_fluid(i) };
Else
wall() += { bnd_fluid(i) };
EndIf

EndFor

// ======================================================
// 5. PHYSICAL GROUPS
// ======================================================

Physical Surface("fluid") = {400};
Physical Surface("solid") = {300};

Physical Line("cylinder") = {bnd_solid()};
Physical Line("inlet1") = {inlet1()};
Physical Line("inlet2") = {inlet2()};
Physical Line("outlet1") = {outlet1()};
Physical Line("outlet2") = {outlet2()};
Physical Line("wall") = {wall()};
Physical Point("A") = {A};
Physical Point("B") = {B};



// ======================================================
// 6. MESH REFINEMENT: gradual coarsening from particle
// ======================================================

Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeExtendFromBoundary = 0;

Mesh.Algorithm = 6;          // Frontal-Delaunay
Mesh.Smoothing = 30;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

// ------------------------------------------------------
// Mesh sizes
// ------------------------------------------------------

// Particle radius = 8e-6 m.
// h_particle = R/3 gives about 19 elements around the circumference
// and about 6 elements across the diameter.
h_particle = cylinder_radius / 2.5;

// Far-field mesh size.
// channel_height/10 gives about 10 elements across the channel height.
h_far = channel_height / 6;

// Smooth transition length.
// 30R is smoother than 20R or 25R, so the mesh does not jump suddenly.
r_grow = 8 * cylinder_radius;

// ------------------------------------------------------
// Distance from particle boundary
// ------------------------------------------------------

Field[1] = Distance;
Field[1].CurvesList = {bnd_solid()};
Field[1].NumPointsPerCurve = 300;

// Gradual transition from h_particle to h_far
Field[2] = Threshold;
Field[2].InField = 1;

Field[2].DistMin = 0;
Field[2].DistMax = r_grow;

Field[2].SizeMin = h_particle;
Field[2].SizeMax = h_far;

// ------------------------------------------------------
// Keep the particle interior fine
// ------------------------------------------------------

Field[3] = Box;
Field[3].VIn  = h_particle;
Field[3].VOut = h_far;

Field[3].XMin = particle_x_offset - 1.2*cylinder_radius;
Field[3].XMax = particle_x_offset + 1.2*cylinder_radius;
Field[3].YMin = particle_y_offset - 1.2*cylinder_radius;
Field[3].YMax = particle_y_offset + 1.2*cylinder_radius;
Field[3].ZMin = -1;
Field[3].ZMax =  1;

// ------------------------------------------------------
// Use the finest requested size
// ------------------------------------------------------

Field[4] = Min;
Field[4].FieldsList = {2,3};

Background Field = 4;
