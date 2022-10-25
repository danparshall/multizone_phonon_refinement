function elastic = estimate_elastic_line(DAT);

WINDOW = 5;

inds_elastic = find( (DAT.eng > -WINDOW) & (DAT.eng < WINDOW) );
assert(length(inds_elastic) > 1, 'Need finite number of energies for elastic window.');
elastic = max(DAT.y_dat(inds_elastic, :));