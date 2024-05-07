%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example code for cross section as in Fig. 4. Must call run.m in the 
% code directory before using.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
load('cmap')

filename = '../code/frames/f-10.mat'; %change number for corresponding timestamp
render(filename,cmap(1,:))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute contour along cross section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function render(source,c)

  load(source,'x','z','phi')
  [xx,zz] = meshgrid(x,z); xx = xx'; zz = zz'; 
  phi = squeeze(phi(:,end/2,:));
  [b,h] = contour(xx,zz,phi,[0.5 0.5]);
  delete(h)
  fill(b(1,2:end),b(2,2:end),c,'LineWidth',2)
  axis equal off
  xlim([1 6])
  ylim([0 1.5])
  xlabel(''),ylabel('')
  xticklabels(''),yticklabels('')
  drawnow
end