function startvars = empirical_starting_heights(DAT, start_cens, start_wids,  reswid, varargin);
% Find approximate starting heights, based on the measured data.  For use when predictions aren't available.
%
% For each peak center, find the measured intensity at the closest energy.  Scale by number of nearby centers,
% to handle overlap.


if ~exist('reswid')
    reswid = 2.5;
end


% for each center, find index of closest energy
eng = DAT.eng(:);
cens = start_cens(:);
A = repmat(eng, [1 length(cens)]);
[min_val, closest_idx] = min(abs(A-cens'));


% for each center, count number of neighbors (we reduce intensity by 1/neighbors)
n_close = [];
for i_cen = 1:length(cens)
    lo_neighbors = (cens > cens(i_cen) - reswid);
    hi_neighbors = (cens < cens(i_cen) + reswid);
    neighbors = sum(lo_neighbors .* hi_neighbors);
    n_close = [n_close neighbors];
end
scaling = 1./n_close;


% for each Q, find intensity at points closest to our peak centers (and scale)
n_q = size(DAT.y_dat, 2);
obs_intensity = DAT.y_dat(closest_idx, :);
start_hts = repmat(scaling(:), 1, n_q) .* obs_intensity;


% replace any bad heights with a reasonable guess
%avg_ht = median(start_hts(~isnan(start_hts)));
avg_ht = median(median(DAT.y_dat(~isnan(DAT.y_dat))));
start_hts(~(start_hts > 0)) = avg_ht;

startvars = [start_cens(:), start_wids(:), start_hts];