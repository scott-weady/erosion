%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes stl files from field data
% For an example, after running a simulation call
%     iso2stl(fin) 
% with fin an output .mat file (e.g. fin = '../code/frames/f-0.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iso2stl(fin)
  fout = 'mudlion.stl';
  load(fin,'phi')
  [F,V] = isosurface(smooth3(phi),0.5);
  stlwrite(fout,F,V)
end