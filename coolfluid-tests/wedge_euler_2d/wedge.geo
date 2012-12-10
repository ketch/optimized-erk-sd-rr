// This file contains all the gmsh commands to generate a grid made with 
// quadrilateral cells, i.e. quads.
////////////////////////////////////////////////////////////////////////////////

// Parameters
l = 1;  // base and hight of the triangles 
charLength = 1; // characteristic grid size on a edge

x_min = -5.0;
x_max = 15.0;
y_min = -5.0;
y_max = 5.0;

z = 0.0;

//////////////////////////////////////

// Points: triangle's verticies
p1t = newp; Point(p1t) = {-0.5,0.0,z,charLength};
p2t = newp; Point(p2t) = {0.5,-0.5,z,charLength};
p3t = newp; Point(p3t) = {0.5,0.5,z,charLength};

// Points: inlet boundary points
p1in = newp; Point(p1in) = {x_min,y_max,z,charLength};
p2in = newp; Point(p2in) = {x_min,0.0,z,charLength};
p3in = newp; Point(p3in) = {x_min,y_min,z,charLength};

// Points: outlet boundary points
p1out = newp; Point(p1out) = {x_max,y_max,z,charLength};
p2out = newp; Point(p2out) = {x_max,2.0,z,charLength};
p3out = newp; Point(p3out) = {x_max,-2.0,z,charLength};
p4out = newp; Point(p4out) = {x_max,y_min,z,charLength};

// Points: top boundary points (this boundary shares p1in with the inlet and 
// p1out with the outlet)
p1top = newp; Point(p1top) = {-0.5,y_max,z,charLength};
p2top = newp; Point(p2top) = {0.5,y_max,z,charLength};

// Points: bottom boundary points (this boundary shares p3in with the inlet and 
// p3out with the outlet)
p1bot = newp; Point(p1bot) = {-0.5,y_min,z,charLength};
p2bot = newp; Point(p2bot) = {0.5,y_min,z,charLength};


// Lines: triangle
l1t = newl; Line(l1t) = {p1t,p2t};
l2t = newl; Line(l2t) = {p2t,p3t};
l3t = newl; Line(l3t) = {p3t,p1t};

// Lines: inlet 
l1in = newl; Line(l1in) = {p1in,p2in};
l2in = newl; Line(l2in) = {p2in,p3in};

// Lines: outlet 
l1out = newl; Line(l1out) = {p4out,p3out};
l2out = newl; Line(l2out) = {p3out,p2out};
l3out = newl; Line(l3out) = {p2out,p1out};

// Lines: top 
l1top = newl; Line(l1top) = {p1out,p2top};
l2top = newl; Line(l2top) = {p2top,p1top};
l3top = newl; Line(l3top) = {p1top,p1in};

// Lines: bottom 
l1bot = newl; Line(l1bot) = {p3in,p1bot};
l2bot = newl; Line(l2bot) = {p1bot,p2bot};
l3bot = newl; Line(l3bot) = {p2bot,p4out};

// Lines: internal connections
l1con = newl; Line(l1con) = {p2in,p1t};
l2con = newl; Line(l2con) = {p1t,p1top};
l3con = newl; Line(l3con) = {p3t,p2top};
l4con = newl; Line(l4con) = {p3t,p2out};
l5con = newl; Line(l5con) = {p2t,p3out};
l6con = newl; Line(l6con) = {p2t,p2bot};
l7con = newl; Line(l7con) = {p1t,p1bot};

// Line loops
ll1 = newl; Line Loop(ll1) = {l3top,l1in,l1con,l2con};
ll2 = newl; Line Loop(ll2) = {l2top,-l2con,-l3t,l3con};
ll3 = newl; Line Loop(ll3) = {l1top,-l3con,l4con,l3out};
ll4 = newl; Line Loop(ll4) = {-l4con,-l2t,l5con,l2out};
ll5 = newl; Line Loop(ll5) = {-l5con,l6con,l3bot,l1out};
ll6 = newl; Line Loop(ll6) = {-l1t,l7con,l2bot,-l6con};
ll7 = newl; Line Loop(ll7) = {-l1con,l2in,l1bot,-l7con};

// Surfaces
s1 = news; Plane Surface(s1) = {ll1};
s2 = news; Plane Surface(s2) = {ll2};
s3 = news; Plane Surface(s3) = {ll3};
s4 = news; Plane Surface(s4) = {ll4};
s5 = news; Plane Surface(s5) = {ll5};
s6 = news; Plane Surface(s6) = {ll6};
s7 = news; Plane Surface(s7) = {ll7};



// Transfinite lines
Transfinite Line{l3top,-l1con,-l1bot} = 40 Using Progression 1.01;
Transfinite Line{l2top,l3t,-l1t,l2bot} = 40 Using Progression 1.03;
Transfinite Line{-l1top,l4con,l5con,l3bot} = 200 Using Progression 1.01;
Transfinite Line{l2t,l2out} = 40 Using Bump 0.40;
Transfinite Line{-l1in,l2con,l3con,l3out,l2in,l7con,l6con,-l1out} = 60 Using Progression 1.03;

// Transfinite surfaces
Transfinite Surface{s1} = {p1in,p2in,p1t,p1top};
Transfinite Surface{s2} = {p1top,p1t,p3t,p2top};
Transfinite Surface{s3} = {p2top,p3t,p2out,p1out};
Transfinite Surface{s4} = {p3t,p2t,p3out,p2out};
Transfinite Surface{s5} = {p2t,p2bot,p4out,p3out};
Transfinite Surface{s6} = {p1t,p1bot,p2bot,p2t};
Transfinite Surface{s7} = {p2in,p3in,p1bot,p1t};

// Recombine surfaces
Recombine Surface(s1);
Recombine Surface(s2);
Recombine Surface(s3);
Recombine Surface(s4);
Recombine Surface(s5);
Recombine Surface(s6);
Recombine Surface(s7);

// Physical lines
Physical Line("triangle") = {3, 1, 2};    // Triangle
Physical Line("inlet") = {4, 5};          // Inlet
Physical Line("outlet") = {8, 7, 6};      // Outlet
Physical Line("top") = {11, 10, 9};       // Top
Physical Line("bottom") = {12, 13, 14};   // Bottom

// Physical surface
Physical Surface("domain") = {29, 30, 31, 32, 33, 34, 35}; // Infield
