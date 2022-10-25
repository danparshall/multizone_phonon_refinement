function [modelout,jac_out,varargout]=calc_singleQ(AUX,idx);
% [modelout,jac_out,varargout]=calc_singleQ(AUX,DAT,idx);
% 	Calculates the data for a single Q.  
%	If requested, calculates and delivers Jacobian as well.
%	Cycles through, calculates profile and jacobian for each phonon.

if ~exist('idx'); idx=1; end

eng = AUX.eng(AUX.mask(:,idx))';     % row vector of just energies being fitted

cen = AUX.auxvars(1:end-1, 1, 1);
hts = AUX.auxvars(1:end-1, idx+1, 1);
wid = AUX.auxvars(1:end-1, 1, 2);
res = AUX.auxvars(1:end-1, 1+idx, 2);


%%%% WIDTHS
FWHM = AUX.auxvars(1:end-1, 1, 2) + AUX.auxvars(1:end-1, idx+1, 2); % add phonon+resolution width
asymm = AUX.peak_asymmetry;
w1 = FWHM * asymm/(asymm+1);
w2 = FWHM * 1/(asymm+1);


% Background
constant = AUX.auxvars(end,idx+1,1);
slope = AUX.auxvars(end,idx+1,2);
%[cen, wid, hts, res]
%[constant, slope]

%%%% INIT
modelout = zeros(size(eng));

if nargout==2;
	% number of columns of jacobian is equal to 4*nrows of freevars.
	% Each batch corresponds to:
	%	cens + 1 dummy
	%	wids + 1 dummy
	%	heights + constant BG
	%	reswids + linear BG
	vsize = size(AUX.freevars, 1);
	jac_out = zeros(length(eng), 4*vsize);
	jac_out(:, 3*vsize) = ones(length(eng),1);	% const BG
	jac_out(:, 4*vsize) = eng(:);				% linear BG
end


%%%% Calculate profile for each phonon
for ind=1:AUX.Nph
	if nargout==2	% when asked for jacobian
		[splitgauss, splitjac] = calc_splitgauss_full(eng, cen(ind), hts(ind), wid(ind), res(ind), asymm );
		ind_cwh = [0, 2*vsize, vsize] + ind;  	% resorts so that first Nph cols are cens, then wids, then heights
		jac_out(:, ind_cwh) = splitjac;

	else			% no jacobian
		splitgauss = calc_splitgauss_full(eng, cen(ind), hts(ind), wid(ind), res(ind), asymm );
	end

	modelout = modelout + splitgauss';
end

% tack on the backgrounds to the convolved splitgausses
modelout = modelout + constant;		% constant BG
modelout = modelout + (slope*eng);	% linear BG
