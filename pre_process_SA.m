% Bryan Herman
% Pre-process script for discontinuity factor example

% create regions here
reg(1).fluxDX(1) = 37.3529; % volume integrated flux g=1
reg(1).fluxDX(2) = 10.7199; % volume integrated flux g=2
reg(1).totRR(1) = 20.6895; % total reaction rate g=1
reg(1).totRR(2) = 14.1948; % total reaction rate g=2
reg(1).scattRR(1,1) = 19.6758; % scattering reaction rate h=1-->g=1 
reg(1).scattRR(2,1) = 1.36679E-02; % scattering reaction rate h=2-->g=1
reg(1).scattRR(1,2) = 0.689451; % scattering reaction rate h=1-->g=2
reg(1).scattRR(2,2) = 13.5054; % scattering reaction rate h=2-->g=2
reg(1).nfissRR(1,1) = 0.176007; % fission reaction rate h=1-->g=1 
reg(1).nfissRR(2,1) = 0.808498; % fission reaction rate h=2-->g=1
reg(1).nfissRR(1,2) = 0.0; % fission reaction rate h=1-->g=2
reg(1).nfissRR(2,2) = 0.0; % fission reaction rate h=2-->g=2
reg(1).current_l(1) = 0.0; % current on left surface g=1
reg(1).current_l(2) = 0.0; % current on left surface g=2
reg(1).current_r(1) = 0.0; % current on right surface g=1
reg(1).current_r(2) = 0.0; % current on right surface g=2
reg(1).diffusion(1) = 181.811; % diffusion coefficient in group 1
reg(1).diffusion(2) = 11.6054; % diffusion coefficient in group 2
reg(1).keff = 0.98451; % multiplication factor

% heterogeneous flux at surface
flux_het_r1 = 1.10933;
flux_het_r2 = 0.297646;

% compute ADFs
ADF_r1 = (flux_het_r1/reg(1).fluxDX(1))*(34/1);
ADF_r2 = (flux_het_r2/reg(1).fluxDX(2))*(34/1);

reg(2).fluxDX(1) = 36.8940; % volume integrated flux g=1
reg(2).fluxDX(2) = 6.71652; % volume integrated flux g=2
reg(2).totRR(1) = 20.3771; % total reaction rate g=1
reg(2).totRR(2) =  9.01911; % total reaction rate g=2
reg(2).scattRR(1,1) = 19.3647; % scattering reaction rate h=1-->g=1 
reg(2).scattRR(2,1) = 1.24062E-02; % scattering reaction rate h=2-->g=1
reg(2).scattRR(1,2) = 0.644479; % scattering reaction rate h=1-->g=2
reg(2).scattRR(2,2) = 8.37463; % scattering reaction rate h=2-->g=2
reg(2).nfissRR(1,1) = 0.265966; % fission reaction rate h=1-->g=1 
reg(2).nfissRR(2,1) = 0.959094; % fission reaction rate h=2-->g=1
reg(2).nfissRR(1,2) = 0.0; % fission reaction rate h=1-->g=2
reg(2).nfissRR(2,2) = 0.0; % fission reaction rate h=2-->g=2
reg(2).current_l(1) = 0.0; % current on left surface g=1
reg(2).current_l(2) = 0.0; % current on left surface g=2
reg(2).current_r(1) = 0.0; % current on right surface g=1
reg(2).current_r(2) = 0.0; % current on right surface g=2
reg(2).diffusion(1) = 179.492; % diffusion coefficient in group 1
reg(2).diffusion(2) = 7.12936; % diffusion coefficient in group 2
reg(2).keff = 1.22506; % multiplication factor

% flux left surface
flux_het_l1 = 1.09164;
flux_het_l2 = 0.182014;

% compute ADFs
ADF_l1 = (flux_het_l1/reg(2).fluxDX(1))*(34/1);
ADF_l2 = (flux_het_l2/reg(2).fluxDX(2))*(34/1);

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
        keff = reg(i).keff;
        
        % compute balance
        bal(i,g) = leakage + inters - scatt - 1/keff*fiss;
        
    end
    
end

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

% perform upscatter corrected downscattering
%reg(1).scattxs(1,2) = reg(1).scattxs(1,2) - reg(1).scattxs(2,1)*reg(1).fluxDX(2)/reg(1).fluxDX(1);
%reg(2).scattxs(1,2) = reg(2).scattxs(1,2) - reg(1).scattxs(2,1)*reg(2).fluxDX(2)/reg(2).fluxDX(1);

% now zero out downscattering xs
% reg(1).scattxs(2,1) = 0;
% reg(2).scattxs(2,1) = 0;