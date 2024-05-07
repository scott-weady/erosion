%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface and erosion rate visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function renderSurface(xx,yy,zz,phi,Vn,shade)

  % Set up lighting
  h1 = camlight; lightangle(h1,-90,40)
  h2 = camlight; lightangle(h2,0,40)
  h3 = camlight; lightangle(h3,90,40)

  % Rearrange
  xx = permute(xx,[2 1 3]);
  yy = permute(yy,[2 1 3]);
  zz = permute(zz,[2 1 3]);
  phi = smooth3(permute(phi,[2 1 3])); %smooth for visualization
  Vn = smooth3(permute(Vn,[2 1 3])); % " "

  % Show object
  if shade == true
    [f,v,c] = isosurface(xx,yy,zz,phi,0.5,Vn);
    patch('Vertices',v,'Faces',f,'FaceVertexCData',c, ...
          'FaceColor','interp','EdgeColor', 'none')
  else
    s = patch(isosurface(xx,yy,zz,phi,0.5));
    s.FaceColor = [0.588 0.294 0];
    s.EdgeColor = 'none';
    s.FaceLighting = 'gouraud';
  end

  % Custom colormap
  load('cmap','cmap')
  colormap(cmap)

  material dull
  axis equal off
  view([-30 20])

end
