% Bryan Herman
% Pre-process script for discontinuity factor example

% create regions here
reg(1).fluxDX(1) = 13.2231; % volume integrated flux g=1
reg(1).fluxDX(2) = 3.83637; % volume integrated flux g=2
reg(1).totRR(1) = 7.37725; % total reaction rate g=1
reg(1).totRR(2) =  5.07601; % total reaction rate g=2
reg(1).scattRR(1,1) = 7.00975; % scattering reaction rate h=1-->g=1 
reg(1).scattRR(2,1) = 4.99069E-03; % scattering reaction rate h=2-->g=1
reg(1).scattRR(1,2) = 0.250804; % scattering reaction rate h=1-->g=2
reg(1).scattRR(2,2) = 4.82963; % scattering reaction rate h=2-->g=2
reg(1).nfissRR(1,1) = 6.24632E-02; % fission reaction rate h=1-->g=1 
reg(1).nfissRR(2,1) = 0.288782; % fission reaction rate h=2-->g=1
reg(1).nfissRR(1,2) = 0.0; % fission reaction rate h=1-->g=2
reg(1).nfissRR(2,2) = 0.0; % fission reaction rate h=2-->g=2
reg(1).current_l(1) = 0.0; % current on left surface g=1
reg(1).current_l(2) = 0.0; % current on left surface g=2
reg(1).current_r(1) = 0.187380 - 0.241610; % current on right surface g=1
reg(1).current_r(2) = 4.85344E-02 - 4.41094E-02; % current on right surface g=2
reg(1).diffusion(1) = 63.7680; % diffusion coefficient in group 1
reg(1).diffusion(2) = 4.15852; % diffusion coefficient in group 2
reg(1).keff = 1.13934; % multiplication factor

reg(2).fluxDX(1) = 23.8350; % volume integrated flux g=1
reg(2).fluxDX(2) = 4.31347; % volume integrated flux g=2
reg(2).totRR(1) = 13.1124; % total reaction rate g=1
reg(2).totRR(2) =  5.79609; % total reaction rate g=2
reg(2).scattRR(1,1) = 12.4670; % scattering reaction rate h=1-->g=1 
reg(2).scattRR(2,1) = 7.88653E-03; % scattering reaction rate h=2-->g=1
reg(2).scattRR(1,2) = 0.410058; % scattering reaction rate h=1-->g=2
reg(2).scattRR(2,2) = 5.38161; % scattering reaction rate h=2-->g=2
reg(2).nfissRR(1,1) = 0.171005; % fission reaction rate h=1-->g=1 
reg(2).nfissRR(2,1) = 0.617085; % fission reaction rate h=2-->g=1
reg(2).nfissRR(1,2) = 0.0; % fission reaction rate h=1-->g=2
reg(2).nfissRR(2,2) = 0.0; % fission reaction rate h=2-->g=2
reg(2).current_l(1) = 0.187380 - 0.241610; % current on left surface g=1
reg(2).current_l(2) = 4.85344E-02 - 4.41094E-02; % current on left surface g=2
reg(2).current_r(1) = 0.0; % current on right surface g=1
reg(2).current_r(2) = 0.0; % current on right surface g=2
reg(2).diffusion(1) = 116.653; % diffusion coefficient in group 1
reg(2).diffusion(2) = 4.57234; % diffusion coefficient in group 2
reg(2).keff = 1.13934; % multiplication factor

% compute balance in each group
nregs = length(reg);
ng = length(reg(1).fluxDX);

% preallocate balance matrix
bal = zeros(nregs,ng);

% begin loop around regions
for i = 1:nregs
    
    % begin loop around groups
    for g = 1:ng
    
        % zero out accumulated vars
        scatt = 0;
        fiss = 0;
        
        % compute net leakage
        leakage = reg(i).current_r(g) - reg(i).current_l(g);
        
        % get interactions
        inters = reg(i).totRR(g);
        
        % loop around inscatter groups
        for h = 1:ng
            
            % append scattering
            scatt = scatt + reg(i).scattRR(h,g);
            fiss = fiss + reg(i).nfissRR(h,g);
            
        end
 
        % get keff
        keff = reg(1).keff;
        
        % compute balance
        bal(i,g) = leakage + inters - scatt - 1/keff*fiss;
        
    end
    
end

% record actual k
k_OPENMC = keff;

% compute cross sections
reg(1).totxs = reg(1).totRR./reg(1).fluxDX;
reg(1).scattxs(:,1) = reg(1).scattRR(:,1)./reg(1).fluxDX';
reg(1).scattxs(:,2) = reg(1).scattRR(:,2)./reg(1).fluxDX';
reg(1).nfissxs(:,1) = reg(1).nfissRR(:,1)./reg(1).fluxDX';
reg(1).nfissxs(:,2) = reg(1).nfissRR(:,2)./reg(1).fluxDX';
reg(1).diff = reg(1).diffusion./reg(1).fluxDX;

reg(2).totxs = reg(2).totRR./reg(2).fluxDX;
reg(2).scattxs(:,1) = reg(2).scattRR(:,1)./reg(2).fluxDX';
reg(2).scattxs(:,2) = reg(2).scattRR(:,2)./reg(2).fluxDX';
reg(2).nfissxs(:,1) = reg(2).nfissRR(:,1)./reg(2).fluxDX';
reg(2).nfissxs(:,2) = reg(2).nfissRR(:,2)./reg(2).fluxDX';
reg(2).diff = reg(2).diffusion./reg(2).fluxDX;