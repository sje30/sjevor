%% Compute Voronoi areas of a set of mosaics, both dmin and other
%% areas.  These functions use matlab_vor.m to do the hard work.

%% Thu 27 Feb 2003
for n = 1:6
  in = sprintf('dmin%d.txt', n)
  out = sprintf('dmin%d_matareas.txt', n)
  matlab_vor(in, out);
end
matlab_vor('w81s.on.d', 'w81s.on.matareas.txt')
matlab_vor('triarray.dat', 'triarray_matareas.txt')

