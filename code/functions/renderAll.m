%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function renderAll(xx,yy,zz,u,v,w,phi,Vn)

  shade = false; %set to true to color by erosion rate
  renderSurface(xx,yy,zz,phi,Vn,shade) %show object
  renderStreamlines(xx,yy,zz,u,v,w) %show streamlines

  % Set visualization domain
  Lx = max(xx(:))+min(xx(:));
  Ly = max(yy(:))+min(yy(:));
  Lz = max(zz(:))+min(zz(:));
  axis equal
  axis([0 Lx 0 Ly 0 Lz])

  % Add background
  hold on
  fill3([0 Lx Lx 0 0],[0 0 Ly Ly 0],[0 0 0 0 0],0.5*[1 1 1],'facealpha',0.5)
  fill3([Lx Lx Lx Lx Lx],[0 0 Ly Ly 0],[0 Lz Lz 0 0],0.5*[1 1 1],'facealpha',0.25)
  fill3([0 Lx Lx 0 0],[Ly Ly Ly Ly Ly],[0 0 Lz Lz 0],0.5*[1 1 1],'facealpha',0.25)

  % Format figure
  fig = gcf;
  fig.Units = 'inches';
  fig.PaperSize = [8 4];
  fig.PaperPosition = [0 0 fig.PaperSize];
  fig.Position = [4 4 0 0] + fig.PaperPosition;
  fig.Color = 'k';
  fig.InvertHardcopy = 'off';
  
  % Format axes
  ax = gca;
  ax.Units = 'inches';
  ax.Position(1:2) = 0;
  ax.Position(3:4) = fig.PaperSize;
  ax.Color = 'k';

  view([-20 5]);

end