% Test
%% Geometry Setup
geom.map = [1 2];
geom.dx = [22.6 12.0];
geom.mesh = [500 500];
geom.totmesh = cumsum(geom.mesh);
geom.albedo = [1,1];
geom.f_l = [[1;1],[1;1]];
geom.f_r = [[1;1],[1;1]];
geom.current_r = [[0.123494893307097;-0.003968964699779],[0.0;0.0]];
geom.current_l = [[0.0;0.0],[0.123494893307097;0.003968964699779]];

%% Material Setup

% UO2 Assembly
mat(1).totxs = [0.312167,1.49284];
mat(1).scattxs = [0.297949,0.0051751,4.89388e-5,0.349585];
mat(1).nfissxs = [0.0221029,1.658,0,0];
mat(1).keff = 1.3662;
mat(1).diff = [1/(3*0.253273),1/(3*1.41404)];

% UO2 Assembly
mat(2).totxs = [0.320625,0.370885];
mat(2).scattxs = [0.313324,0.00218517,5.2433e-04,0.345478];
mat(2).nfissxs = [0.0019385,0.0166938,0,0];
mat(2).keff = 1.3662;
mat(2).diff = [1/(3*0.269713),1/(3*0.330466)];

%% Run Problem

% Create Loss Matrix
M = setup_M(geom,mat,'albedo');

% Create Production Matrix
F = setup_F(geom,mat);

% Solve for eigenvalue/eigenvector
[k_both,phi_both] = solve_problem(M,F);

% Compute Discontinuity Factors
[geom,phimat] = compute_discontinuity(geom,mat);

% Create Loss Matrix
M = setup_M(geom,mat,'albedo');

% Create Production Matrix
F = setup_F(geom,mat);

% Solve for eigenvalue/eigenvector
[k_both_DF,phi_both_DF] = solve_problem(M,F);