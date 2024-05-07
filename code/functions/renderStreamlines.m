%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Streamline visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function renderStreamlines(xx,yy,zz,u,v,w)
  
  % Get domain size
  Lx = max(xx(:))+min(xx(:));
  Ly = max(yy(:))+min(yy(:));
  Lz = max(zz(:))+min(zz(:));

  Ns = 10; %number of streamlines
  h = xx(2)-xx(1); %grid spacing
  sx = 2*h; %x-location of initial seed
  sy = 0.95*Ly/2; %y-location of initial seed
  sz = linspace(0,Lz,Ns); sz = sz(2:end-1); %z-location of initial seed

  [sx,sy,sz] = meshgrid(sx,sy,sz); %create streamline grid
  sl = streamline(stream3(permute(xx,[2 1 3]),...
      permute(yy,[2 1 3]),...
      permute(zz,[2 1 3]),...
      permute(u,[2 1 3]),...
      permute(v,[2 1 3]),...
      permute(w,[2 1 3]),...
      sx,sy,sz));

  set(sl,'LineWidth',3,'Color',[1 0.655 0.188])

end