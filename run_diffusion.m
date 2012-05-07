% Bryan Herman
% Simple 1-D Diffusion Code

%% Reset
clear
close all

%% Compute balance for both
pre_process_both

%% Geometry Setup
geom.map = [1 2];
geom.dx = [21.42 21.42];
geom.mesh = [50 50];
geom.totmesh = cumsum(geom.mesh);
geom.albedo = [1,1];
geom.f_l = [[1;1],[1;1]];
geom.f_r = [[1;1],[1;1]];
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
% mat(1).diff = [reg(1).diff(1),reg(1).diff(2)];
mat(1).keff = reg(1).keff;
mat(1).diff = [1/(3*reg(1).totxs(1)),1/(3*reg(1).totxs(2))];
mat(1).diff(1) = 1.1;
mat(1).diff(2) = 0.3;

% UO2 Assembly
mat(2).totxs = [reg(2).totxs(1),reg(2).totxs(2)];
mat(2).scattxs = [reg(2).scattxs(1,1),reg(2).scattxs(1,2), ...
                  reg(2).scattxs(2,1),reg(2).scattxs(2,2)];
mat(2).nfissxs = [reg(2).nfissxs(1,1),reg(2).nfissxs(1,2), ...
                  reg(2).nfissxs(2,1),reg(2).nfissxs(2,2)];
% mat(2).diff = [reg(2).diff(1),reg(2).diff(2)];
mat(2).keff = reg(2).keff;
mat(2).diff = [1/(3*reg(2).totxs(1)),1/(3*reg(2).totxs(2))];
mat(2).diff(1) = 1.1;
mat(2).diff(2) = 0.3;

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

%% Compute balance for SA
pre_process_SA

%% Geometry Setup
geom.map = [1 2];
geom.dx = [21.42 21.42];
geom.mesh = [50 50];
geom.totmesh = cumsum(geom.mesh);
geom.albedo = [1,1];
geom.f_l = [[1;1],[1;1]];
geom.f_r = [[1;1],[1;1]];
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
% mat(1).diff = [reg(1).diff(1),reg(1).diff(2)];
mat(1).keff = reg(1).keff;
mat(1).diff = [1/(3*reg(1).totxs(1)),1/(3*reg(1).totxs(2))];
mat(1).diff(1) = 1.1;
mat(1).diff(2) = 0.3;

% UO2 Assembly
mat(2).totxs = [reg(2).totxs(1),reg(2).totxs(2)];
mat(2).scattxs = [reg(2).scattxs(1,1),reg(2).scattxs(1,2), ...
                  reg(2).scattxs(2,1),reg(2).scattxs(2,2)];
mat(2).nfissxs = [reg(2).nfissxs(1,1),reg(2).nfissxs(1,2), ...
                  reg(2).nfissxs(2,1),reg(2).nfissxs(2,2)];
% mat(2).diff = [reg(2).diff(1),reg(2).diff(2)];
mat(2).keff = reg(2).keff;
mat(2).diff = [1/(3*reg(2).totxs(1)),1/(3*reg(2).totxs(2))];
mat(2).diff(1) = 1.1;
mat(2).diff(2) = 0.3;
%% Run Problem

% Create Loss Matrix
M = setup_M(geom,mat,'albedo');

% Create Production Matrix
F = setup_F(geom,mat);

% Solve for eigenvalue/eigenvector
[k_SA,phi_SA] = solve_problem(M,F);


%% Rerun SA with DFs
geom.f_l = [[1;1],[ADF_l1;ADF_l2]];
geom.f_r = [[ADF_r1;ADF_r2],[1;1]];

% Create Loss Matrix
M = setup_M(geom,mat,'albedo');

% Create Production Matrix
F = setup_F(geom,mat);

% Solve for eigenvalue/eigenvector
[k_SA_DFs,phi_SA_DFs] = solve_problem(M,F);

%% Plot all
plot_all(phimat,phi_SA,phi_both,phi_both_DF,phi_SA_DFs);
