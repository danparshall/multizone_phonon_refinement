function refinedCenters = remove_degenerate(rawCenters)

refinedCenters = [];
while length(rawCenters) != 0
	refinedCenters = [refinedCenters rawCenters(1)];
	mask = find(rawCenters != rawCenters(1));
	rawCenters = rawCenters(mask);
end
