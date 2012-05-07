function [geom,phimat] = compute_discontinuity(geom,mat)

% get length of dimensions
nx = length(geom.map); % spatial mesh
ng = length(mat(1).totxs); % energy mesh

for imesh = 1:nx
    
   % create an aritificial geometry object with only 1 coarse mesh
   geomnew.map = geom.map(imesh);
   geomnew.dx = geom.dx(imesh);
   geomnew.mesh = geom.mesh(imesh);
   geomnew.totmesh = cumsum(geomnew.mesh);
   geomnew.f_l = [[1;1],[1;1]];
   geomnew.f_r = [[1;1],[1;1]];
   geomnew.current_l = geom.current_l(:,imesh);
   geomnew.current_r = geom.current_r(:,imesh);

   % compute mesh width
   dx = geomnew.dx/geomnew.mesh;
   
   % compute loss matrix
   M = setup_M(geomnew,mat,'current');
   
   % compute production matrix
   F = setup_F(geomnew,mat);
   
   % create righthandside
   Q = zeros(size(F,1),1);
   
   % begin loop around groups
   for g = 1:ng
       
       % get row index for left cell
       row_l = indices_to_matrix(g,1,ng);
       
       % set left group current
       Q(row_l) = geomnew.current_l(g);
       
       % get row index for right cell
       row_r = indices_to_matrix(g,geomnew.mesh,ng);
       
       % set right group current
       Q(row_r) = -geomnew.current_r(g);
       
   end
   
   % divide Q by dx
   Q = Q/dx;
   
   % compute homogeneous flux
   phi = (M - 1/mat(geomnew.map).keff*F)\Q;

   % begin loop around groups
   for g = 1:ng
      
       % get row indices
       row_l = indices_to_matrix(g,1,ng);
       row_r = indices_to_matrix(g,geomnew.mesh,ng);
       
       % get diffusion coefficient for cell
       diff = mat(geomnew.map).diff(g);
       
       % compute left homogeneous surface flux
       phi_ls = phi(row_l) + (geomnew.current_l(g)*dx)/(2*diff);
       
       % compute right homogeneous surface flux
       phi_rs = phi(row_r) - (geomnew.current_r(g)*dx)/(2*diff);
       
       % compute left discontinuity factor
       geom.f_l(g,imesh) = 1.0/phi_ls;
       
       % compute right discontinuity factor
       geom.f_r(g,imesh) = 1.0/phi_rs;
       
   end
   
   % bank flux in matrix
   phimat{imesh} = reshape(phi,ng,geomnew.mesh);
   
end

% renormalize DFs
for imesh = 1:nx-1
    
    for g = 1:ng
        
        rat = geom.f_r(g,imesh)/geom.f_l(g,imesh+1);
        geom.f_l(g,imesh+1) = 2/(1+rat);
        geom.f_r(g,imesh) = rat*geom.f_l(g,imesh+1);
        
    end
    
end

end


function matidx = indices_to_matrix(g,i,ng)

    % compute matrix location
    matidx = g + ng*(i-1);
    
end
