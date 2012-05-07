function M = setup_M(geom,mat,BC)

% compute indices
nx = geom.totmesh(length(geom.totmesh)); % spatial mesh
ng = length(mat(1).totxs); % energy mesh
n = nx*ng; % total dimensions of matrix
nnz = (nx-2)*3*ng + 2*2*ng + (ng^2-ng)*nx; % number of nonzeros

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

    % get cell data
    diff_c = mat(matidx).diff(g);
    dx = geom.dx(meshidx)/geom.mesh(meshidx);

    % set all DFs to 1
    f_rite_l = 1;
    f_left_l = 1;
    f_rite_r = 1;
    f_left_r = 1;
    
    % check for left boundary
    if i == 1
        
        % compute dtilde for boundary cell left
        if strcmp(BC,'albedo')
            alb_l = geom.albedo(1);
            dtilde_l = (2*diff_c*(1-alb_l))/(4*diff_c*(1+alb_l) + (1-alb_l)*dx);
        elseif strcmp(BC,'current')
            dtilde_l = 0;
        else
            error('No BC specified!')
        end
        
    else
        
        % get diff of left cell
        meshidx_l = find(i-1<=geom.totmesh,1,'first');
        matidx_l = geom.map(meshidx_l);
        diff_l = mat(matidx_l).diff(g);
        dx_l = geom.dx(meshidx_l)/geom.mesh(meshidx_l);
        
        % get discontinuity factors around interface
        if meshidx ~= meshidx_l
            f_rite_l = geom.f_r(g,meshidx_l); % left cell right DF
            f_left_l = geom.f_l(g,meshidx);   % right cell left DF
        end
        
        % compute dtilde for interior cell
        dtilde_l = (2*diff_c*diff_l)/(f_rite_l*dx_l*diff_c + f_left_l*dx*diff_l);
        
        % get matrix column and bank value
        icol = indices_to_matrix(g,i-1,ng);
        rows(kount) = irow;
        cols(kount) = icol;
        vals(kount) = -f_rite_l*dtilde_l/dx;
        kount = kount + 1;
        
    end
    
    % check for right boundary
    if i == nx
        
        % compute dtilde for boundary cell right
        if strcmp(BC,'albedo')
            alb_r = geom.albedo(2);
            dtilde_r = (2*diff_c*(1-alb_r))/(4*diff_c*(1+alb_r) + (1-alb_r)*dx);
        elseif strcmp(BC,'current')
            dtilde_r = 0;
        end
        
    else
        
        % get diff of right cell
        meshidx_r = find(i+1<=geom.totmesh,1,'first');
        matidx_r = geom.map(meshidx_r);
        diff_r = mat(matidx_r).diff(g);
        dx_r = geom.dx(meshidx_r)/geom.mesh(meshidx_r);
        
        % get discontinuity factors around interface
        if meshidx ~= meshidx_r
            f_rite_r = geom.f_r(g,meshidx);   % left cell right DF
            f_left_r = geom.f_l(g,meshidx_r); % right cell left DF
        end
        
        % compute dtilde for interior cell
        dtilde_r = (2*diff_c*diff_r)/(f_left_r*dx_r*diff_c + f_rite_r*dx*diff_r);
        
        % get matrix column and bank value
        icol = indices_to_matrix(g,i+1,ng);
        rows(kount) = irow;
        cols(kount) = icol;
        vals(kount) = -f_left_r*dtilde_r/dx;
        kount = kount + 1;
        
    end
    
    % calculate diagonal term
    gidx = g + ng*(g - 1);
    remxs = mat(matidx).totxs(g) - mat(matidx).scattxs(gidx);
    leakin = f_rite_r*dtilde_r/dx + f_left_l*dtilde_l/dx;
    rows(kount) = irow;
    cols(kount) = irow;
    vals(kount) = leakin + remxs;
    kount = kount + 1;
    
    % begin loop for scattering
    for h = 1:ng
       
        % cycle if g==h
        if g==h
            continue
        end
        
        % get column number
        icol = indices_to_matrix(h,i,ng);
        
        % get xs
        hidx = g + ng*(h-1);
        scattxshg = mat(matidx).scattxs(hidx);
        rows(kount) = irow;
        cols(kount) = icol;
        vals(kount) = -scattxshg;
        kount = kount + 1;
        
    end
     
end

% create sparse matrix
M = sparse(rows,cols,vals);

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