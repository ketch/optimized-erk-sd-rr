// This file contains all the gmsh commands to generate a grid made with 
// quadrilateral cells, i.e. quads.
////////////////////////////////////////////////////////////////////////////////

// Parameters
l = 1;  // base and hight of the triangles 
nb_cells = 10;
unit = l/nb_cells;
x_min = -5.0;
x_max = 15.0;
y_min = -5.0;
y_max = 5.0;

bl_h  = 1.0;
w_l  = 1.5;

z = 0.0;

//////////////////////////////////////

// Points: triangle's verticies
p1t = newp; Point(p1t) = {-0.5,0.0,z,unit};
p2t = newp; Point(p2t) = {0.5,0.5,z,unit};
p3t = newp; Point(p3t) = {0.5,-0.5,z,unit};


// Points boundary layer
// top left
p1bl = newp; Point(p1bl) = {-0.5,bl_h,z,unit};
// top right
p2bl = newp; Point(p2bl) = {0.5,0.5+bl_h,z,unit};
// bottom left
p3bl = newp; Point(p3bl) = {-0.5,-bl_h,z,unit};
// bottom right
p4bl = newp; Point(p4bl) = {0.5,-0.5-bl_h,z,unit};


// Points attack
p1a = newp; Point(p1a) = {-0.5-bl_h/2,bl_h,z,unit};
p2a = newp; Point(p2a) = {-0.5-bl_h/2,0.0,z,unit};
p3a = newp; Point(p3a) = {-0.5-bl_h/2,-bl_h,z,unit};


// Points wake
// top
p1w = newp; Point(p1w) = {0.5+w_l,0.5+w_l,z,2*unit};
// top middle
p2w = newp; Point(p2w) = {0.5+w_l,0.5,z,unit};
// bottom middle
p3w = newp; Point(p3w) = {0.5+w_l,-0.5,z,unit};
// bottom
p4w = newp; Point(p4w) = {0.5+w_l,-0.5-w_l,z,2*unit};


// Points: inlet boundary points
// top
p1in = newp; Point(p1in) = {x_min,y_max,z,10*unit};
// top middle
p2in = newp; Point(p2in) = {x_min,bl_h,z,10*unit};
// middle
p3in = newp; Point(p3in) = {x_min,0,z,10*unit};
// bottom middle
p4in = newp; Point(p4in) = {x_min,-bl_h,z,10*unit};
// bottom
p5in = newp; Point(p5in) = {x_min,y_min,z,10*unit};

// Points: outlet boundary points
p1out = newp; Point(p1out) = {x_max,y_max,z,5*unit};
p2out = newp; Point(p2out) = {x_max,0.5+w_l,z,5*unit};
p3out = newp; Point(p3out) = {x_max,0.5,z,5*unit};
p4out = newp; Point(p4out) = {x_max,-0.5,z,5*unit};
p5out = newp; Point(p5out) = {x_max,-0.5-w_l,z,5*unit};
p6out = newp; Point(p6out) = {x_max,y_min,z,5*unit};


// Points: top boundary points (this boundary shares p1in with the inlet and 
// p1out with the outlet)
p1top = p1in;
p2top = newp; Point(p2top) = {-0.5,y_max,z,10*unit};
p3top = newp; Point(p3top) = {0.5,y_max,z,10*unit};
p4top = newp; Point(p4top) = {0.5+bl_h,y_max,z,10*unit};
p5top = p1out;

// Points: bottom boundary points (this boundary shares p3in with the inlet and 
// p3out with the outlet)
p1bot = p5in;
p2bot = newp; Point(p2bot) = {-0.5,y_min,z,10*unit};
p3bot = newp; Point(p3bot) = {0.5,y_min,z,10*unit};
p4bot = newp; Point(p4bot) = {0.5+bl_h,y_min,z,10*unit};
p5bot = p6out;


// Lines: top boundary
l1top = newl; Line(l1top) = {p1top,p2top};
l2top = newl; Line(l2top) = {p2top,p3top};
l3top = newl; Line(l3top) = {p3top,p4top};
l4top = newl; Line(l4top) = {p4top,p5top};

// Lines: inlet boundary
l1in = newl; Line(l1in) = {p1in,p2in};
l2in = newl; Line(l2in) = {p2in,p3in};
l3in = newl; Line(l3in) = {p3in,p4in};
l4in = newl; Line(l4in) = {p4in,p5in};

// Lines: bottom boundary
l1bot = newl; Line(l1bot) = {p1bot,p2bot};
l2bot = newl; Line(l2bot) = {p2bot,p3bot};
l3bot = newl; Line(l3bot) = {p3bot,p4bot};
l4bot = newl; Line(l4bot) = {p4bot,p5bot};

// Line: outlet boundary
l1out = newl; Line(l1out) = {p1out,p2out};
l2out = newl; Line(l2out) = {p2out,p3out};
l3out = newl; Line(l3out) = {p3out,p4out};
l4out = newl; Line(l4out) = {p4out,p5out};
l5out = newl; Line(l5out) = {p5out,p6out};

// Lines: connection between boundary, BL, and triangle
l1hor1  = newl; Line(l1hor1)  = {p1a,p1bl};
l2hor1  = newl; Line(l2hor1)  = {p1bl,p2bl};
l3hor1  = newl; Line(l3hor1)  = {p2bl,p1w};
l4hor1  = newl; Line(l4hor1)  = {p1w,p2out};
                          
l1hor2  = newl; Line(l1hor2)  = {p2a,p1t};
l2hor2  = newl; Line(l2hor2)  = {p1t,p2t};
l3hor2  = newl; Line(l3hor2)  = {p2t,p2w};
l4hor2  = newl; Line(l4hor2)  = {p2w,p3out};

l1hor3  = l1hor2;
l2hor3  = newl; Line(l2hor3)  = {p1t,p3t};
l3hor3  = newl; Line(l3hor3)  = {p3t,p3w};
l4hor3  = newl; Line(l4hor3)  = {p3w,p4out};

l1hor4  = newl; Line(l1hor4)  = {p3a,p3bl};
l2hor4  = newl; Line(l2hor4)  = {p3bl,p4bl};
l3hor4  = newl; Line(l3hor4)  = {p4bl,p4w};
l4hor4  = newl; Line(l4hor4)  = {p4w,p5out};


l1a  = newl; Line(l1a)  = {p1a,p2a};
l2a  = newl; Line(l2a)  = {p2a,p3a};

//l1ver1  = newl; Line(l1ver1)  = {p2top,p1bl};
l2ver1  = newl; Line(l2ver1)  = {p1bl,p1t};
l3ver1  = newl; Line(l3ver1)  = {p1t,p3bl};
//l4ver1  = newl; Line(l4ver1)  = {p3bl,p2bot};

//l1ver2  = newl; Line(l1ver2)  = {p3top,p2bl};
l2ver2  = newl; Line(l2ver2)  = {p2bl,p2t};
l3ver2  = newl; Line(l3ver2)  = {p2t,p3t};
l4ver2  = newl; Line(l4ver2)  = {p3t,p4bl};
//l5ver2  = newl; Line(l5ver2)  = {p4bl,p3bot};

//l1ver3  = newl; Line(l1ver3)  = {p4top,p1w};
l2ver3  = newl; Line(l2ver3)  = {p1w,p2w};
l3ver3  = newl; Line(l3ver3)  = {p2w,p3w};
l4ver3  = newl; Line(l4ver3)  = {p3w,p4w};
//l5ver3  = newl; Line(l5ver3)  = {p4w,p4bot};

// Line loops


ll14 = newl; Line Loop(ll14) = {l2hor2,-l2ver2,-l2hor1,l2ver1};
ll15 = newl; Line Loop(ll15) = {l3hor2,-l2ver3,-l3hor1,l2ver2};
ll16 = newl; Line Loop(ll16) = {l3hor3,-l3ver3,-l3hor2,l3ver2};

ll17 = newl; Line Loop(ll17) = {l3hor4,-l4ver3,-l3hor3,l4ver2};
ll18 = newl; Line Loop(ll18) = {l2hor4,-l4ver2,-l2hor3,l3ver1};

ll1a = newl; Line Loop(ll1a) = {l1hor2,-l2ver1,-l1hor1,l1a};
ll2a = newl; Line Loop(ll2a) = {l1hor4,-l3ver1,-l1hor3,l2a};

// Line loops
// top left to top right
llround = newl; Line Loop(llround) = {l1bot,l2bot,l3bot,l4bot,-l5out,-l4hor4,-l3hor4,-l2hor4,-l1hor4,-l2a,-l1a,l1hor1,l2hor1,l3hor1,l4hor1,-l1out,-l4top,-l3top,-l2top,-l1top,l1in,l2in,l3in,l4in};
ll1a = newl; Line Loop(ll1a) = {l1hor2,-l2ver1,-l1hor1,l1a};
ll2a = newl; Line Loop(ll2a) = {l1hor4,-l3ver1,-l1hor3,l2a};
ll1w  = newl; Line Loop(ll1w) = {l4hor2,-l2out,-l4hor1,l2ver3};
ll2w  = newl; Line Loop(ll2w) = {l4hor3,-l3out,-l4hor2,l3ver3};
ll3w  = newl; Line Loop(ll3w) = {l4hor4,-l4out,-l4hor3,l4ver3};

//llw = newl; Line Loop(llw) = {l4hor4,-l4out,-l3out,-l2out,-l4hor1,l2ver3,l3ver3,l4ver3};

// Surfaces
sround = news; Plane Surface(sround) = {llround};
s1w = news; Plane Surface(s1w) = {ll1w};
s2w = news; Plane Surface(s2w) = {ll2w};
s3w = news; Plane Surface(s3w) = {ll3w};
//sw = news; Plane Surface(sw) = {llw};
s1a = news; Plane Surface(s1a) = {ll1a};
s2a = news; Plane Surface(s2a) = {ll2a};

// Surfaces
s14 = news; Plane Surface(s14) = {ll14};
s15 = news; Plane Surface(s15) = {ll15};
s16 = news; Plane Surface(s16) = {ll16};
s17 = news; Plane Surface(s17) = {ll17};
s18 = news; Plane Surface(s18) = {ll18};


// Transfinite lines
cells = (nb_cells)/2+1;
Transfinite Line{l1hor1,l1hor2,l1hor3,l1hor4} = cells Using Progression 1;
Transfinite Line{-l2hor2,-l2hor3} = cells Using Bump 0.2;
Transfinite Line{-l2hor1,-l2hor4} = cells;

Transfinite Line{l3hor2,l3hor3} = 2*cells Using Progression 1.1;
Transfinite Line{l3hor1,l3hor4} = 2*cells;

Transfinite Line{l3ver2,l3ver3} = 2*cells Using Bump 0.25;
Transfinite Line{l3ver3} = 2*cells;
Transfinite Line{-l1a,l2a,-l2ver3,l4ver3} = 2*cells Using Progression 1.1;
Transfinite Line{-l2ver1,-l2ver2, l3ver1,l4ver2} = 2*cells Using Progression 1.1;

Transfinite Line{l4hor1,l4hor2,l4hor3,l4hor4} = 5*cells Using Progression 1.05;
Transfinite Line{-l2out,l4out} = 2*cells Using Progression 1.1;
Transfinite Line{l3out} = 2*cells;

// Structured part of the mesh (mainly BLs)
Transfinite Surface{s1a} = {p2a,p1t,p1bl,p1a};
Transfinite Surface{s2a} = {p3a,p3bl,p1t,p2a};

Transfinite Surface{s1w} = {p2w,p3out,p2out,p1w};
Transfinite Surface{s2w} = {p3w,p4out,p3out,p2w};
Transfinite Surface{s3w} = {p4w,p5out,p4out,p3w};

Transfinite Surface{s14} = {p1t,p2t,p2bl,p1bl};
Transfinite Surface{s15} = {p2t,p2w,p1w,p2bl};
Transfinite Surface{s16} = {p3t,p3w,p2w,p2t};
Transfinite Surface{s17} = {p4bl,p4w,p3w,p3t};
Transfinite Surface{s18} = {p3bl,p4bl,p3t,p1t};


// Recombine surfaces

Recombine Surface(sround);
Recombine Surface(s1w);
Recombine Surface(s2w);
Recombine Surface(s3w);
//Recombine Surface(sw);
Recombine Surface(s1a);
Recombine Surface(s2a);
Recombine Surface(s14);
Recombine Surface(s15);
Recombine Surface(s16);
Recombine Surface(s17);
Recombine Surface(s18);


// Physical lines (BC)
Physical Line("inlet")    = {l1in, l2in, l3in, l4in};      // Inlet
Physical Line("top")      = {l1top, l2top, l3top, l4top};      // Top boundary
Physical Line("outlet")   = {l1out, l2out, l3out, l4out, l5out};  // Outlet
Physical Line("bottom")   = {l1bot, l2bot, l3bot, l4bot};      // Bottom boundary
Physical Line("triangle") = {l2hor2, l3ver2, l2hor3};             // Triangle

// Physical surface (computational domain)
//Physical Surface("domain") = {sround,sw,s1a,s2a,s14,s15,s16,s17,s18};
Physical Surface("domain") = {sround,s1w,s2w,s3w,s1a,s2a,s14,s15,s16,s17,s18};

// Algorithm type (pave in gambit)
Mesh.SubdivisionAlgorithm=1;
