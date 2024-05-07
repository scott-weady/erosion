
function fig = journal_figure(sz,sc,tex,fontname)

if nargin < 2
  sc = 1;
end
if nargin < 3
  tex = 'latex';
end

% set(gcf,'defaultAxesFontSize',12*sc)
set(groot,'defaultlinelinewidth',1.25*sc)
set(groot,'defaultAxesLineWidth',1*sc)
set(groot,'defaultAxesFontSize',8*sc)
set(groot,'defaulttextinterpreter',tex); 
set(groot,'defaultAxesTickLabelInterpreter',tex); 
set(groot,'defaultLegendInterpreter',tex);
if strcmp(tex,'tex')
  set(groot,'defaultaxesFontName',fontname)
end

fig = figure;
fig.Units = 'inches';
fig.PaperSize = sz*sc;
fig.PaperPosition = [0 0 sz]*sc;
fig.Position = [2 2 sz*sc];
fig.PaperPositionMode = 'auto';
