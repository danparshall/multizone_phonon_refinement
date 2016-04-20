function [mask, goodheight, center_height] = make_mask(SYMDAT,startvars,sqwind)
% [mask, goodheight, free_iCenht] = make_mask(SYMDAT,startvars,sqwind)
%	looks through data, excludes datasets that don't have a good peak
%	mask is a boolean, size(ydat). 1 for valid data.
%	goodheight is finite for good peaks, 0 otherwise
%	free_iCenht is (Nph x Nq+1), indicates which phonons/heights may be fit


%User editable values
peak_expansion = 1;% multiple of resolution FWHM out from peak center to consider "peak"
BG_expansion = 3;% multiple of resolution FWHM out from peak center to consider as background (excludes peak)
peak_fraction = .4;% of all the possible places y values can exist near a peak, the minimum fraction that must exist to let that peak pass
BG_fraction = .1;% of all the possible places y values can exist around (but not on) a peak, the minimum fraction that must exist to let that peak pass
peak2error = 2;% ratio of peak height to error bar where if the error bar is to big, the peak is not counted...
						%(set peak2error to any negative number to disable exclusion of proportionally high-error data)

%starting variables
SQW = SYMDAT{sqwind}.SQW;
centers = startvars(:,1);
eng = SQW.xdat(:,1);
xstep = eng(2)-eng(1);
ydat = SQW.ydat;
edat = SQW.edat;

Nph=length(centers);
Nq=size(ydat,2);


%excludes all y values where there is no data
mask=~isnan(ydat);

%finds index of energy level closest to peak center
eindex = [];
for ind = 1:Nph
	thisval=find(abs(eng-centers(ind))<=xstep/2);
	if isempty(thisval)
		disp(' WARNING in "make_mask" : phonon center out of energy range');
	else
		eindex = [eindex thisval];
	end
end


% eBins is a list of indices that are the boundary point of the bins around each phonon 
eBins = [1];
for ind = 1:length(eindex)-1
	halfwayPoint = mean([eindex(ind) eindex(ind+1)]);
	eBins = [eBins halfwayPoint halfwayPoint];
end
eBins = [eBins length(eng)];
eBins = round(eBins);
peaksize = [];%peakSize is the distance (in meV) from a peak that the mask looks (should be within a magnitude of resolution width probably)
for ind = 1:Nph
	peaksize = [peaksize peak_expansion*merchop(SYMDAT{sqwind}.Ei,SYMDAT{sqwind}.chopfreq,centers(ind))];
end


%removes peaks that aren't dense enough or the error is too big and creates height fitting matrices
ebar = ydat./edat;
ebar(find(abs(ebar)<peak2error)) = NaN;
ebar(isnan(ebar)) = 0;

goodheight = zeros(Nph,Nq);
heights_free = ones(Nph,Nq);
centers_free = ones(Nph,1);

for iCen = 1:Nph
	%creates index for area on peak
	if peaksize(iCen)<eng(eindex(iCen))-eng(eBins(1))
		edexL1 = eindex(iCen)-round(peaksize(iCen)/xstep);
	else
		edexL1 = eBins(iCen*2-1);
	end
	if peaksize(iCen)<eng(eindex(iCen))-eng(eBins(end))
		edexH1 = eindex(iCen)+round(peaksize(iCen)/xstep);
	else
		edexH1 = eBins(iCen*2);
	end
	%creates index for area around, but not on, peak
	if 2*peaksize(iCen)<eng(eindex(iCen))-eng(eBins(1))
		edexL2 = eindex(iCen)-round(peaksize(iCen)*2/xstep);
	else
		edexL2 = eBins(iCen*2-1);
	end
	if 2*peaksize(iCen)<eng(eindex(iCen))-eng(eBins(end))
		edexH2 = eindex(iCen)+round(peaksize(iCen)*2/xstep);
	else
		edexH2 = eBins(iCen*2);
	end

	peak_possible = length(edexL1:edexH1);% maximum possible number of points on the "peak"
	BG_possible = length([edexL2:edexL1,edexH1:edexL2]);% maximum possible number of points in the "background"

	for pho = 1:Nq
		%excludes any data that has too low a density of points on a peak/ near a peak
		peak_points = sum(mask(edexL1:edexH1,pho));% actual # of good points on peak
		BGpoints = sum(mask([edexL2:edexL1,edexH1:edexL2],pho));% actual # of good points in BG
		if  peak_points/peak_possible < peak_fraction | BGpoints/BG_possible < BG_fraction% are there not enough actual points?
			mask(eBins(iCen*2-1):eBins(iCen*2),pho) = 0;
			heights_free(iCen,pho) = 0;
		else
			goodheight(iCen,pho) = max(ydat(eBins(iCen*2-1):eBins(iCen*2),pho));
		end
	end

	%excludes data where the errorbars are too big relative to the peak size
	mask(eBins(iCen*2-1):eBins(iCen*2),find(sum(ebar(edexL1:edexH1,:))==0)) = 0;
	heights_free(iCen,find(sum(ebar(edexL1:edexH1,:))==0)) = 0;

	%peaks with no data at all should be ignored
	if sum(sum(mask(eBins(iCen*2-1):eBins(iCen*2),:))) == 0
		centers_free(iCen) = 0;
	end
end
%this is where it is decided whether a center/height is to be fit
goodheight = goodheight.*heights_free;%heights that don't get fit are set to zero
center_height = [centers_free heights_free];
