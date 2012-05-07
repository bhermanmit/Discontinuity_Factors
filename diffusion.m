% Bryan Herman
% Simple 1-D Diffusion Code

%% Geometry Setup
geom.map = [1 2];
geom.dx = [20 20];
geom.mesh = [20 20];
geom.totmesh = cumsum(geom.mesh);
geom.albedo = [1,1];
geom.fleft = [[1;1],[1;1]];
geom.fright = [[1;1],[1;1]];
geom.current_r = [[reg(1).current_r(1);reg(1).current_r(2)] ...
                 [reg(2).current_r(1);reg(2).current_r(2)]];
geom.current_l = [[reg(1).current_l(1);reg(1).current_l(2)] ...
                 [reg(2).current_l(1);reg(2).current_l(2)]];

%% Material Setup

% UO2 Assembly
mat(1).totxs = [reg(1).totxs(1),reg(1).totxs(2)];
mat(1).scattxs = [reg(1).scattxs(1,1),reg(1).scattxs(1,2), ...
                  reg(1).scattxs(2,1),reg(1).scattxs(2,2)];
mat(1).nfissxs = [reg(1).nfissxs(1,1),reg(1).nfissxs(1,2), ...
                  reg(1).nfissxs(2,1),reg(1).nfissxs(2,2)];
mat(1).diff = [reg(1).diff(1),reg(1).diff(2)];

% UO2 Assembly
mat(2).totxs = [reg(2).totxs(1),reg(2).totxs(2)];
mat(2).scattxs = [reg(2).scattxs(1,1),reg(2).scattxs(1,2), ...
                  reg(2).scattxs(2,1),reg(2).scattxs(2,2)];
mat(2).nfissxs = [reg(2).nfissxs(1,1),reg(2).nfissxs(1,2), ...
                  reg(2).nfissxs(2,1),reg(2).nfissxs(2,2)];
mat(2).diff = [reg(2).diff(1),reg(2).diff(2)];

%% Run Problem

% Compute Discontinuity Factors
geom = compute_discontinuity(geom,mat);

% Create Loss Matrix
M = setup_M(geom,mat,'albedo');

% Create Production Matrix
F = setup_F(geom,mat);

% Solve for eigenvalue/eigenvector
[k,phi] = solve_problem(M,F);