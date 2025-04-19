// Cube dimensions
L = 1.0;  // Length of the cube

// Points
Point(1) = {0, 0, 0, L/10};
Point(2) = {L, 0, 0, L/10};
Point(3) = {L, L, 0, L/10};
Point(4) = {0, L, 0, L/10};
Point(5) = {0, 0, L, L/10};
Point(6) = {L, 0, L, L/10};
Point(7) = {L, L, L, L/10};
Point(8) = {0, L, L, L/10};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {1, 5};
Line(6) = {2, 6};
Line(7) = {3, 7};
Line(8) = {4, 8};
Line(9) = {5, 6};
Line(10) = {6, 7};
Line(11) = {7, 8};
Line(12) = {8, 5};

// Surfaces
Line Loop(1) = {1, 2, 3, 4};  // Bottom face
Plane Surface(1) = {1};
Line Loop(2) = {5, 9, -6, -1};  // Front face
Plane Surface(2) = {2};
Line Loop(3) = {6, 10, -7, -2};  // Right face
Plane Surface(3) = {3};
Line Loop(4) = {7, 11, -8, -3};  // Back face
Plane Surface(4) = {4};
Line Loop(5) = {8, 12, -5, -4};  // Left face
Plane Surface(5) = {5};
Line Loop(6) = {9, 10, 11, 12};  // Top face
Plane Surface(6) = {6};

// Volume
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

// Physical groups
Physical Volume("conductor") = {1};
Physical Surface("dirichlet_0") = {1};    // z=0 surface (anode)
Physical Surface("dirichlet_1") = {6};    // z=L surface (cathode)
Physical Surface("neumann_0") = {2, 3, 4, 5};  // side surfaces (insulator)

// Mesh settings
Mesh.Algorithm = 6;  // Frontal-Delaunay for 3D
Mesh.Algorithm3D = 1;  // Delaunay
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

// Save settings
Mesh.SaveAll = 1;  // Save all elements (including physical groups)
Mesh.Binary = 0;   // Save in ASCII format
Mesh.SaveNodes = 1;  // Save nodes
Mesh.SaveElements = 1;  // Save elements
Mesh.SaveGroupsOfNodes = 1;  // Save node groups
Mesh.SaveGroupsOfElements = 1;  // Save element groups
Mesh.SaveElementTagType = 1;  // Save element tag type 