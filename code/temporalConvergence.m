%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporal convergence study for phase-field erosion model 
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

Nz = 64; %grid resolution
tf = 1; %final time
tsave = inf; %save time (off)
tplt = inf; %plotting time (off)
dt = 0.02; %initial time step

% Begin refinement
main; phim1 = phi;
for nrefine = 1:Nrefine
  dt = dt/2; main;
  dphi = abs(phi(:)-phim1(:));
  l1 = mean(dphi); %l1 norm
  linf = max(dphi); %linf norm
  err(nrefine,:) = [dt,l1,linf]; %save errors
  phim1 = phi; %prepare for next refinement
end

%% Plot errors
close all
loglog(err(:,1),err(:,2),'o-','DisplayName','$||\phi_{\Delta t} - \phi_{\Delta t/2}||_1$'), hold on
loglog(err(:,1),err(:,3),'o-','DisplayName','$||\phi_{\Delta t} - \phi_{\Delta t/2}||_\infty$')
loglog(err(:,1),0.5*err(end,2)*err(:,1)/err(end,1),'k--','DisplayName','$\sim \Delta t$')
xlabel('time step, $\Delta t$'), ylabel('refinement error')
legend('location','northwest','edgecolor','none','color','none')
xlim([err(end,1)/2 2*err(1,1)]), ylim([min(min(err(:,2:3)))/4 2*max(max(err(:,2:3)))])
