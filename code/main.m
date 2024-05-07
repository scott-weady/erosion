%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main script for phase-field erosion simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup
close all
addpath('functions')
target = 'frames'; %directory for saving images
if exist(target,'dir'); rmdir(target,'s'); end
mkdir(target)

%% Discretization
Nx = round(Lx/Lz)*Nz; %streamwise resolution
Ny = round(Ly/Lz)*Nz; %horizontal resolution
x = avg(linspace(0,Lx,Nx+1)); hx = Lx/Nx; %x-grid
y = avg(linspace(0,Ly,Ny+1)); hy = Ly/Ny; %y-grid
z = avg(linspace(0,Lz,Nz+1)); hz = Lz/Nz; %z-grid
[yy,xx,zz] = meshgrid(y,x,z); %3d grid
dV = hx*hy*hz; %volume differential

del = 2*hx; %interface damping length 
eta = 0.01*Re*del^2; %penalization parameter
m = 0.2*hx^2; %curvature smoothing coefficient

Nt = floor(tf/dt+0.1); %number of time steps
nplt = floor(tplt/dt+0.1); %plotting interval
nsave = floor(tsave/dt+0.1); %saving interval
nf = 0; %start frame counter
nt = 0; %start time step counter

%% Initial condition
fprintf("Setting up initial conditions...\n")

% Boundary conditions (N - top, S - bottom, E - outlet, W - inlet)
uN = 1; uS = 0; uE = 0; uW = tanh(0.5*zz(1,:,:)/hx);
vN = 0; vS = 0; vE = 0; vW = 0;
wN = 0; wS = 0; wE = 0; wW = 0;

% Preallocate variables
u = uW.*ones(Nx-1,Ny,Nz); %x-velocity
v = zeros(Nx,Ny-1,Nz); %y-velocity
w = zeros(Nx,Ny,Nz-1); %z-velocity

% Principal axes of ellipsoid
lx = 2; %x-coordinate
lz = 3/2; %y-coordinate
ly = 2/3; %z-coordinate

% Center of ellipsoid
x0 = 0.4*Lx; %x (toward
y0 = 0.5*Ly; %y (center of domain)
z0 = -0.01*Lz; %z (shift to avoid vertical wall)

% Create mask function for clay
rr = sqrt((xx - x0).^2/lx^2 + (yy - y0).^2/ly^2 + (zz - z0).^2/lz^2);
phi1 = 0.5*(1 - tanh(2*(rr-1)/del)); %horizontal axis
phi2 = 0.5*(1 - tanh(2*(zz-3)/del)); %vertical axis
phi = phi1.*phi2;

% Create inclusion
r_inc = 1/3; %inclusion radius
h_inc = 1/2; %inclusion height
x_inc = -5/6; %inclusion location
z_inc = 3/4; %inclusion location
rr = sqrt((xx - (x0-(lx-r_inc)-x_inc)).^2 + (yy - y0).^2);

C1 = 0.5*(1 - tanh(2*(rr-r_inc)/del));
C2 = 0.5*(1 - tanh(-2*(zz-z_inc)/del));
C3 = 0.5*(1 - tanh(2*(zz-z_inc-h_inc)/del));
C = C_clay + (C_inc-C_clay)*C1.*C2.*C3;

u = u.*(1-avgx(phi)); %zero velocity inside body
v = v.*(1-avgy(phi)); % " "
w = w.*(1-avgz(phi)); % " "

V0 = sum(phi(:))*dV; V = V0; %initial volume
Vol = [0 V]; %store initial volume
Vnmax = [0 0]; %store max normal velocity

%% Pressure Poisson solver
fprintf("Building Poisson solver...\n")
Lp_h = kron(speye(Ny),L_fd(Nx,hx,1,3)) + kron(L_h(Ny,hy),speye(Nx));
Lp_h = kron(speye(Nz),Lp_h) + kron(L_fd(Nz,hz,1,1),speye(Nx*Ny));
[Rp,~,Pp] = chol(Lp_h); %compute Cholesky decomposition
Lp = struct('R',Rp,'P',Pp,'Rt',Rp.','Pt',Pp.');

%% Main loop
fprintf("Starting main loop...\n")
while nt < Nt && V/V0 > 1e-3
  
  t1 = tic;

  % Add ghost cells
  uu = extend_u(u,uN,uS,uE,uW);
  vv = extend_v(v,vN,vS,vE,vW);
  ww = extend_w(w,wN,wS,wE,wW);
  pphi = extend_phi(phi);

  % Compute rhs of phi equation
  tau = shear(uu,vv,ww,pphi,hx,hy,hz); %tangential velocity |grad(phi) x u|
  fphi = -C.*tau/del + m*(laplacian(pphi,hx,hy,hz) - phi.*(1 - phi).*(1 - 2*phi)/del^2); 

  % Advection (conservative form)
  ua = avgy(uu(:,:,2:end-1));
  va = avgx(vv(:,:,2:end-1));
  uvx = Dx(ua.*va,hx); uvx = uvx(:,2:end-1,:);
  uvy = Dy(ua.*va,hy); uvy = uvy(2:end-1,:,:);
  
  ua = avgz(uu(:,2:end-1,:));
  wa = avgx(ww(:,2:end-1,:));
  uwx = Dx(ua.*wa,hx); uwx = uwx(:,:,2:end-1);
  uwz = Dz(ua.*wa,hz); uwz = uwz(2:end-1,:,:);
  
  va = avgz(vv(2:end-1,:,:));
  wa = avgy(ww(2:end-1,:,:));
  vwy = Dy(va.*wa,hy); vwy = vwy(:,:,2:end-1);
  vwz = Dz(va.*wa,hz); vwz = vwz(:,2:end-1,:);
  
  ua = avgx(uu(:,2:end-1,2:end-1));
  va = avgy(vv(2:end-1,:,2:end-1));
  wa = avgz(ww(2:end-1,2:end-1,:));

  uux = Dx(ua.^2,hx);
  vvy = Dy(va.^2,hy);
  wwz = Dz(wa.^2,hz);

  % Compute rhs of velocity equation
  fu = -(uux + uvy + uwz) + (1/Re)*laplacian(uu,hx,hy,hz) - (1/eta)*avgx(phi).*u;
  fv = -(uvx + vvy + vwz) + (1/Re)*laplacian(vv,hx,hy,hz) - (1/eta)*avgy(phi).*v;
  fw = -(uwx + vwy + wwz) + (1/Re)*laplacian(ww,hx,hy,hz) - (1/eta)*avgz(phi).*w;

  % Advance
  u = u + dt*fu;
  v = v + dt*fv;
  w = w + dt*fw;
  phi = phi + dt*fphi;
  phi = max(min(phi,1),0);

  % Pressure projection
  uu = extend_u(u,uN,uS,uE,uW);
  vv = extend_v(v,vN,vS,vE,vW);
  ww = extend_w(w,wN,wS,wE,wW);
  p = pressure(uu,vv,ww,hx,hy,hz,Lp);
  u = u - Dx(p,hx);
  v = v - Dy(p,hy);
  w = w - Dz(p,hz);
   
  % Compute volume and store
  V = sum(phi(:))*dV; %compute volume
  Vol = [Vol;dt*nt,V]; %store volume

  % Compute erosion rate and store
  Vn = 4*C.*tau; %erosion rate
  Vnmax = [Vnmax;dt*nt,max(Vn(:))]; %save max normal velocity

  % Print out information
  clc; el = toc(t1);
  fprintf('%6s: %1.4f\n%6s: %1.4f\n%6s: %1.4es\n','time',dt*nt,'volume',V/V0,'clock',el)
  if isnan(u(1)); fprintf(' unstable, terminating...\n'); break; end

  % Save volume and max erosion rate
  if mod(nt,100) == 0
    save(strcat(target,'/','V'),'Vol','Vnmax')
  end

  % Visualization
  if mod(nt,nplt) == 0 && nplt < inf

    cla
    renderAll(xx,yy,zz,ua,va,wa,phi,Vn)
    drawnow

    filename = strcat(target,'/',sprintf('f-%d\n',nf));
    if mod(nt,nsave) == 0
      save(filename,'x','y','z','phi','ua','va','wa','tau','C')
    end
    print(filename,'-dpng','-r500')
    nf = nf + 1;

  end
   
  nt = nt + 1;
   
end

%% Functions

% Cell center-face averaging operators
function ua = avg(u)
  ua = 0.5*(u(2:end) + u(1:end-1));
end
function ua = avgx(u)
  ua = 0.5*(u(2:end,:,:) + u(1:end-1,:,:)); 
end
function ua = avgy(u)
  ua = 0.5*(u(:,2:end,:) + u(:,1:end-1,:)); 
end
function uavg = avgz(u)
  uavg = 0.5*(u(:,:,2:end) + u(:,:,1:end-1)); 
end

% First derivative operators
function ux = Dx(u,hx)
  ux = diff(u,[],1)/hx; 
end
function uy = Dy(u,hy)
  uy = diff(u,[],2)/hy; 
end
function uz = Dz(u,hz)
  uz = diff(u,[],3)/hz; 
end

% Laplace operator
function Lu = laplacian(uu,hx,hy,hz)
  uxx = diff(uu,2,1)/hx^2; uxx = uxx(:,2:end-1,2:end-1);
  uyy = diff(uu,2,2)/hy^2; uyy = uyy(2:end-1,:,2:end-1);  
  uzz = diff(uu,2,3)/hz^2; uzz = uzz(2:end-1,2:end-1,:);  
  Lu = uxx + uyy + uzz;
end

% Second-order negative Laplacian
% bc1,bc2: Neumann=1, Dirichlet (face)=2, Dirichlet (center)=3;
function Delta = L_fd(N,h,bc1,bc2)
  Delta = spdiags([-1 bc1 0;ones(N-2,1)*[-1 2 -1];0 bc2 -1],-1:1,N,N)'/h^2;
end
% Diagonalized second-order negative Laplacian (periodic)
function Delta = L_h(N,h)
  Delta = spdiags((2/h)*sin(pi/(N)*[0:(N/2) (-N/2+1):-1]'),0,N,N).^2;
end

% Ghost cell extension operators
function uu = extend_u(u,uN,uS,uE,uW)
  uu = zeros(size(u)+2);
  i = 2:(size(uu,1)-1);
  j = 2:(size(uu,2)-1);
  k = 2:(size(uu,3)-1);
  uu(i,j,k) = u;
  uu(1,j,k) = uW;
  uu(end,j,k) = (4*u(end,:,:) - u(end-1,:,:))/3;
  uu(i,[1 end],k) = u(:,[end 1],:);
  uu(i,j,1) = 2*uS-u(:,:,1);
  uu(i,j,end) = 2*uN-u(:,:,end);
end
function vv = extend_v(v,vN,vS,vE,vW)
  vv = zeros(size(v)+2);
  i = 2:(size(vv,1)-1);
  j = 2:(size(vv,2)-1);
  k = 2:(size(vv,3)-1);
  vv(i,j,k) = v;
  vv(1,j,k) = 2*vW-v(1,:,:);
  vv(end,j,k) = v(end,:,:);
  vv(i,[1 end],k) = v(:,[end 1],:);
  vv(i,j,1) = 2*vS-v(:,:,1);
  vv(i,j,end) = 2*vN-v(:,:,end);
end
function ww = extend_w(w,wN,wS,wE,wW)
  ww = zeros(size(w)+2);
  i = 2:(size(ww,1)-1);
  j = 2:(size(ww,2)-1);
  k = 2:(size(ww,3)-1);
  ww(i,j,k) = w;
  ww(1,j,k) = 2*wW-w(1,:,:);
  ww(end,j,k) = w(end,:,:);
  ww(i,[1 end],k) = w(:,[end 1],:);
  ww(i,j,1) = wS;
  ww(i,j,end) = wN;
end
function pphi = extend_phi(phi)
  pphi = zeros(size(phi)+2);
  i = 2:(size(pphi,1)-1);
  j = 2:(size(pphi,2)-1);
  k = 2:(size(pphi,3)-1);
  pphi(i,j,k) = phi;
  pphi([1 end],j,k) = phi([1 end],:,:);
  pphi(i,[1 end],k) = phi(:,[end 1],:);
  pphi(i,j,[1 end]) = phi(:,:,[1 end]);
end

% Pressure poisson solver
function p = pressure(uu,vv,ww,hx,hy,hz,Lp)
  
  % Pressure projection
  rhs = Dx(uu(:,2:end-1,2:end-1),hx) + ...
        Dy(vv(2:end-1,:,2:end-1),hy) + ...
        Dz(ww(2:end-1,2:end-1,:),hz); %divergence
  shp = size(rhs); %get shape
  rhs_h = fft(rhs,[],2); rhs_h = reshape(rhs_h,[],1); %take fft
  p_h = -Lp.P*(Lp.R\(Lp.Rt\(Lp.Pt*rhs_h))); %solve
  p_h = reshape(p_h,shp); %reshape
  p = real(ifft(p_h,[],2)); %take ifft

end

% Phase-field approximation of shear-stress
function tau = shear(uu,vv,ww,pphi,hx,hy,hz)

  % Centered difference of phi
  phix = avgx(Dx(pphi(:,2:end-1,2:end-1),hx));
  phiy = avgy(Dy(pphi(2:end-1,:,2:end-1),hy));
  phiz = avgz(Dz(pphi(2:end-1,2:end-1,:),hz));

  % Cell-centered values
  ua = avgx(uu(:,2:end-1,2:end-1));
  va = avgy(vv(2:end-1,:,2:end-1));
  wa = avgz(ww(2:end-1,2:end-1,:));

  % Shear stress
  tau_x =  (phiy.*wa - phiz.*va);
  tau_y = -(phix.*wa - phiz.*ua);
  tau_z =  (phix.*va - phiy.*ua);
  tau = sqrt(tau_x.^2 + tau_y.^2 + tau_z.^2);

end