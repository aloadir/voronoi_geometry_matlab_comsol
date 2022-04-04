% This script is based on a project on the Mathworks website* that creates
% Voronoi diagrams. Here I developed a way to use these diagrams from
% Matlab and export the to Comsol, via the Livelink for Matlab module of
% Comsol.

% * https://www.mathworks.com/matlabcentral/fileexchange/50772-polytope-bounded-voronoi-diagram-in-2d-and-3d

% Download the required files to the Voronoi diagram creation from its Github repository
% If you want to automatic download this files just uncoment the following
% lines:
% urlwrite("https://github.com/hyongju/Polytope-bounded-Voronoi-diagram/blob/master/MY_con2vert.m","MY_con2vert.m");
% urlwrite("https://github.com/hyongju/Polytope-bounded-Voronoi-diagram/blob/master/MY_intersect.m","MY_intersect.m");
% urlwrite("https://github.com/hyongju/Polytope-bounded-Voronoi-diagram/blob/master/MY_setdiff.m","MY_setdiff.m");
% urlwrite("https://github.com/hyongju/Polytope-bounded-Voronoi-diagram/blob/master/demo.m","demo.m");
% urlwrite("https://github.com/hyongju/Polytope-bounded-Voronoi-diagram/blob/master/inhull.m","inhull.m");
% urlwrite("https://github.com/hyongju/Polytope-bounded-Voronoi-diagram/blob/master/pbisec.m","pbisec.m");
% urlwrite("https://github.com/hyongju/Polytope-bounded-Voronoi-diagram/blob/master/polybnd_voronoi.m","polybnd_voronoi.m");
% urlwrite("https://github.com/hyongju/Polytope-bounded-Voronoi-diagram/blob/master/vert2lcon.m","vert2lcon.m");

% The original script ("demo.m") that generates the voronoi geometry can create 2 and 
% 3 dimension objectsand the variable "dimension" defines which one. At this
% moment, this code only uses the 2D case. Part of the 3D functions still
% here for future use, if I return to that.
dimension = 2;                  % dimension of the space

% The original script ("demo.m") uses a variable named "m" to define
% a number of points that would be the vertices of the polygon that
% would cointain the voronoi geometry. 

% The original script contained the following lines to generate the random
% boundary points for the voronoi geometry and positions for the cells
% m = 50;                       % number of boundary point-candidates
% bnd0 = rand(m,d);             % generate random boundary point-candidates

% This program substitutes this "m" sided polygon for a square. The size of 
% the sides is defined by the variable "s". Here the program uses a 
% unitary size, and later it will reescale this multiplying it by a scale 
% factor.
s = 1;                          % Size os the square that delimits the geometry 

% The variable "nPoints" defines the number of cells that the voronoi
% geometry will contain
nPoints = 30;                  % number of points. 

tolerance = 1e-15;               % tolerance value used in "inhull.m" (larger value high precision, possible numerical error)

% The variable "pos0" defines random positions for the center of each 
% voronoi cell
pos0 = rand(nPoints,dimension);	% generate random points

% This is a legacy from the original script of the possibility of 2D and 3D
% geometries. Here we are set to use the "case 2", but the "case 3" ooption
% still here for future users.
switch dimension
    case 2
        bnd0 = [0 0; s 0; s s; 0 s];    % square boundaries
    case 3
        bnd0 = rand(m,dimension);       % generate boundary point-candidates
end

% This step is necessary for the more complex polynomial, on the case of
% the 2D square this will only add the [0,0] on the end of the variable
% "bnd0"
K = convhull(bnd0);             % computes the 2-D or 3-D convex hull of the points in matrix bnd0
bnd_pnts = bnd0(K,:);           % take boundary points from vertices of convex polytope formed with the boundary point-candidates

% Check if the random generated points are inside a convex hull in n 
% dimensions. For the 2D case all the points will be valid and this will
% create a array containing an array of ones with lenght of "nPoints"

in = inhull(pos0,bnd0,[],tolerance); 

% We use this function to choose points that are inside the defined boundary
% For the case of the 2D square of size "s" every point will be on the
% boundaries, so the condition "if in(i) == 1 will be true for all points

u1 = 0;
pos = zeros(size(pos0,1),dimension);

for i = 1:size(pos0,1)
    if in(i) == 1
        u1 = u1 + 1;
        pos(u1,:) = pos0(i,:);
    end
end
%% Here the program applies the script that creates the voronoi geometry
% =========================================================================
% INPUTS:
% pos       points that are in the boundary      n x d matrix (n: number of points d: dimension) 
% bnd_pnts  points that defines the boundary     m x d matrix (m: number of vertices for the convex polytope
% boundary d: dimension)
% -------------------------------------------------------------------------
% OUTPUTS:
% vornb     Voronoi neighbors for each generator point:     n x 1 cells
% vorvx     Voronoi vertices for each generator point:      n x 1 cells
% =========================================================================
[vornb,vorvx] = polybnd_voronoi(pos,bnd_pnts);

clearvars -except vorvx pos s nPoints
%% Initializing the Comsol objects
% The code starts creating the first objects os a Comsol simulation, the model and component nodes.

% Create Model Object
model = ModelUtil.create("Model");
% Create 2D Component
comp1 = model.component.create("comp1",true);

%General properties
% The "s_Matlab" parameter is the same of the "s" variable on the Matlab
% code, that defines the size of the square which is the sample generated.
model.param().set("s_Matlab", "1", "Size os the square that delimits the geometry on Matlab");

% Quantity of random directions that the grains of the ceramic will have. 
quant_dir = 10;     
%% Create Voronoi Geometry
g1 = comp1.geom.create("g1",2);
g1.lengthUnit("mm");

square = model.component("comp1").geom("g1").create("sq1", "Square");
square.set("size", s);

% This loop 

for i=1:size(vorvx,2)
    coord_x = []; 
    coord_y = [];
    for j=1:size(vorvx{1,i},1)
        coord_x = [coord_x vorvx{1,i}(j,1)];
        coord_y = [coord_y vorvx{1,i}(j,2)];
    end
    coord = [];
    coord = [coord_x ; coord_y];
    a = g1.feature.create("c"+i,"BezierPolygon");
    a.set("type","open");
    a.set("degree",1);
    a.set("p",coord);
        
    csol = g1.feature.create("csol"+i,"ConvertToSolid").set("keep", false);
    csol.set("repairtoltype", "absolute");
    csol.set("absrepairtol", 1.0E-15);
    csol.selection("input").set("c"+i);
end

% clearvars a inputset j coord coord_x coord_y nPoints
g1.run;