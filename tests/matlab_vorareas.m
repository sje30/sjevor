%% Compute Voronoi areas of a set of mosaics.
%% Thu 27 Feb 2003
for n = 1:6
  in = sprintf('dmin%d.txt', n)
  out = sprintf('dmin%d_matareas.txt', n)
  dat = load(in);
  [v,c] = voronoin(dat);
  
  npts = length(dat);
  areas = zeros(npts,1);
  for i = 1:npts 
    as = v( c{i}, :);
    areas(i) = polyarea( as(:,1), as(:,2));
  end
  save(out, 'areas', '-ascii')
end

%% save -ascii  dopa_sampleb_matlabarea.txt areas 
