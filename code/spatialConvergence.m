%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatial convergence study for phase-field erosion model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nrefine = 4; %number of refinements
err = zeros(Nrefine,3); %preallocate

Re = 100; %Reynolds number
sigma = 100; %erodibility ratio

% Domain size
Lx = 8; %length
Ly = 4; %width
Lz = 4; %height

C_clay = 0.01; %clay erodibility
C_inc = C_clay/sigma; %inclusion erodibility

Nz = 8; %initial grid resolution
tf = 1; %final time
tsave = inf; %save time (off)
tplt = inf; %plotting time (off)
dt = 0.005; %time step

% Begin refinement
main; phim1 = phi;
for nrefine = 1:Nrefine
  Nz = 2*Nz; main;
  phia = (phi(1:2:end,1:2:end,1:2:end)+phi(2:2:end,2:2:end,2:2:end))/2;
  dphi = abs(phia-phim1); dphi = dphi(:);
  l1 = mean(dphi); %l1 norm
  linf = max(dphi); %linf norm
  err(nrefine,:) = [1/Nz,l1,linf]; %save errors
  phim1 = phi; %prepare for next refinement
end

%% Plot errors
close all
loglog(err(:,1),err(:,2),'o-','DisplayName','$||\phi_{h} - \phi_{h/2}||_1$'), hold on
loglog(err(:,1),err(:,3),'s-','DisplayName','$||\phi_{h} - \phi_{h/2}||_\infty$')
loglog(err(:,1),0.5*err(end,2)*(err(:,1)/err(end,1)),'k--','DisplayName','$\sim h$')
xlabel('grid spacing, $h$'), ylabel('refinement error')
legend('location','northwest','edgecolor','none','color','none')
xlim([err(end,1)/2 2*err(1,1)]), ylim([min(min(err(:,2:3)))/4 1])
