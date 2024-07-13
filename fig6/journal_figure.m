
function fig = journal_figure(sz,sc)

if nargin < 2
  sc = 1;
end

set(groot,'defaultlinelinewidth',1.25*sc)
set(groot,'defaultAxesLineWidth',1*sc)
set(groot,'defaultAxesFontSize',8*sc)
set(groot,'defaulttextinterpreter','latex'); 
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');

fig = figure;
fig.Units = 'inches';
fig.PaperSize = sz*sc;
fig.PaperPosition = [0 0 sz]*sc;
fig.Position = [2 2 sz*sc];
fig.PaperPositionMode = 'auto';
