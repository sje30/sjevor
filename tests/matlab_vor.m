function [] = matlab_vor(in, out)
% IN is name of input data file.
% OUT is name of output data file for writing Voronoi areas.
dat = load(in);
[v,c] = voronoin(dat);

npts = length(dat);
areas = zeros(npts,1);
for i = 1:npts 
  as = v( c{i}, :);
  areas(i) = polyarea( as(:,1), as(:,2));
end
save(out, 'areas', '-ascii')


