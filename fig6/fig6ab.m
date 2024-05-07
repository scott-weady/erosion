%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example code for Fig. 6a,b. Must call run.m in the code directory 
% before using.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
fig = journal_figure([3.25 3],2); %setup figure
filename = '../code/frames/f-10.mat'; %change number for corresponding timestamp
load(filename) %load data

cmax = 0.05; %max range of colorbar
dye = [1 0.6 0]; %streamline color

Vn = 4*tau.*C; %compute erosion rate
[xx,yy,zz] = meshgrid(x,y,z); %make grid
ua = permute(ua,[2 1 3]); %rearrange for streamline function
va = permute(va,[2 1 3]); % " "
wa = permute(wa,[2 1 3]); % " "

%%%% Side view %%%
sp1 = subplot(2,1,1);

renderErosionRate(xx,yy,zz,phi,Vn,true);
clim([0 cmax])
view(0,0)

% Create streamlines
sx = 0.1; sy = max(yy(:))/2-0.05; sz = linspace(0,2.5,14);
[sx,sy,sz] = meshgrid(sx,sy,sz);
sl = streamtube(stream3(xx,yy,zz,ua,va,wa,sx,sy,sz),0.025);
set(sl,'EdgeColor',dye,'FaceAlpha',0.5,'EdgeAlpha',0.5)
xlabel(''), ylabel('')

%%%% Top view %%%
sp2 = subplot(2,1,2);

renderErosionRate(xx,yy,zz,phi,Vn,true);
clim([0 cmax])
view(0,90)

% Create streamlines
sx = 0.1; sy = linspace(0,4,16); sz = 0.33;
[sx,sy,sz] = meshgrid(sx,sy,sz);
sl = streamtube(stream3(xx,yy,zz,ua,va,wa,sx,sy,sz),0.025);
set(sl,'EdgeColor',dye,'FaceAlpha',0.5,'EdgeAlpha',0.5)
xlabel(''),ylabel('')

%%%% Format %%%
load('data/cmap')
colormap(cmap)
clim([0 cmax])

sp1.Units = 'inches';
sp2.Units = 'inches';
sp1.Color = 'k';
sp2.Color = 'k';
sp1.LineWidth = 2;
sp2.LineWidth = 2;

sp2.Position(2) = 0.775;
sp2.Position(3) = 0.85*fig.PaperSize(1);
sp2.Position(4) = sp2.Position(3)*diff(ylim)/diff(xlim);
sp2.Position(1) = (fig.PaperSize(1)-sp2.Position(3))/2;
sp1.Position(2) = sum(sp2.Position([2 4]))+0.2;
sp1.Position([1 3 4]) = sp2.Position([1 3 4]);

hb = colorbar;
hb.Location = 'southoutside';
hb.Units = 'inches';
hb.FontSize = 16;
hb.TickLabelInterpreter = 'latex';
ylabel(hb,'erosion rate, $v_n/U$','interpreter','latex','FontSize',16)
hb.Ticks = 0:0.01:0.1;
hb.Position(2) = 0.6;
hb.Position([1 3]) = sp1.Position([1 3])+0.01*[-1 2];
hb.Position(4) = 0.125;

fig.Color = 'w';
set(fig,'InvertHardCopy','off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create isosurface and color by erosion rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function renderErosionRate(xx,yy,zz,phi,vn,shade)

  % Rearrange and smooth for plotting
  phi = permute(phi,[2 1 3]);
  vn = permute(vn,[2 1 3]);
  phi = smooth3(phi);
  vn = smooth3(vn);
  clay = [0.3922 0.1961 0]; %color of solid

  % Render
  if shade == true
    [f,v,c] = isosurface(xx,yy,zz,phi,0.5,vn);
    patch('Vertices',v,'Faces',f,'FaceVertexCData',c, ...
      'FaceColor','interp','EdgeColor', 'none')
  else
    s = patch(isosurface(xx,yy,zz,phi,0.5));
    s.FaceColor = clay;
    s.EdgeColor = 'none';
    s.FaceLighting = 'gouraud';
  end

  % Format
  material dull
  axis equal
  box on
  xticks([]); yticks([]); zticks([])
  Lx = max(xx(:)); Ly = max(yy(:)); H = 1.75; %plotting range
  axis([Lx/2-2.5 Lx/2+1.5 Ly/2-H/2 Ly/2+H/2 0 H])

end