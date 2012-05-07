function [k,phi] = solve_problem(M,F)

% run eigs function
[phivecs,keigs] = eigs(F,M);

% grab fundamental mode
k = keigs(1,1);
phi = phivecs(:,1);

end