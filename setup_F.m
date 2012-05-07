function F = setup_F(geom,mat)

% compute indices
nx = geom.totmesh(length(geom.totmesh)); % spatial mesh
ng = length(mat(1).totxs); % energy mesh
n = nx*ng; % total dimensions of matrix
nnz = ng^2*nx; % number of nonzeros

% prealloacate arrays
rows = zeros(nnz,1);
cols = zeros(nnz,1);
vals = zeros(nnz,1);
kount = 1;

% begin loop around rows
for irow = 1:n
    
    % get group and x indicies
    [g,i] = matrix_to_indices(irow,nx,ng);
    
    % look up coarse center mesh for materials
    meshidx = find(i<=geom.totmesh,1,'first');
    
    % look up material in mesh idx
    matidx = geom.map(meshidx);
    
    % begin loop for scattering
    for h = 1:ng
        
        % get column number
        icol = indices_to_matrix(h,i,ng);
        
        % get xs
        hidx = g + ng*(h-1);
        nfissxshg = mat(matidx).nfissxs(hidx);
        rows(kount) = irow;
        cols(kount) = icol;
        vals(kount) = nfissxshg;
        kount = kount + 1;
        
    end
     
end

% create sparse matrix
F = sparse(rows,cols,vals);

end

function [g,i] = matrix_to_indices(irow,nx,ng)

    % compute indices
    g = mod(irow-1,ng) +1;
    i = floor(mod(irow-1,ng*nx)/ng) + 1;

end

function matidx = indices_to_matrix(g,i,ng)

    % compute matrix location
    matidx = g + ng*(i-1);
    
end