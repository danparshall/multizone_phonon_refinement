function [modelout,jacout,varargout]=calc_singleQ(AUX,idx);
% [modelout,jacout,varargout]=calc_singleQ(AUX,DAT,idx);
% 	Calculates the data for a single Q.  
%	If requested, calculates and delivers Jacobian as well.
%	Cycles through, calculates profile and jacobian for each phonon.

if ~exist('idx'); idx=1; end

eng = AUX.eng(AUX.mask(:,idx))(:)';     % row vector of just energies being fitted

cen=AUX.auxvars(1:end-1, 1, 1);
ht=AUX.auxvars(1:end-1, idx+1, 1);


%%%% WIDTHS
FWHM = AUX.auxvars(1:end-1, 1, 2) + AUX.auxvars(1:end-1, idx+1, 2); % add phonon+resolution width
asymm = AUX.peak_asymmetry;
w1 = FWHM * asymm/(asymm+1);
w2 = FWHM * 1/(asymm+1);

% Background
constant = AUX.auxvars(end,idx+1,1);
slope = AUX.auxvars(end,idx+1,2);

%%% INIT
modelout = zeros(size(eng));

if nargout==2;
	% number of columns of jacobian is equal to 4*nrows of freevars.
	% Each batch corresponds to:
	%	cens + 1 dummy
	%	wids + 1 dummy
	%	heights + constant BG
	%	reswids + linear BG
	vsize = size(AUX.freevars, 1);
	jacout = zeros(length(eng), 4*vsize);
	jacout(:, 3*vsize) = ones(length(eng),1);	% const BG
	jacout(:, 4*vsize) = eng(:);				% linear BG

%	jacout=[zeros(length(eng),3*AUX.Nph) ones(length(eng),1) eng(:)];  % last two cols are derivatives of background
end


for ind=1:AUX.Nph
	if nargout==2	% when asked for jacobian
		[splitgauss,splitjac] = calc_splitgauss_JAC_fast(eng,cen(ind),ht(ind),w1(ind),w2(ind));
		modelout = modelout + splitgauss';
        
%		ind_j = (3*ind)-2;
%		jacout(:,[ind_j:ind_j+2])=splitjac;     	% cen ht wid (both widths combined)
		ind_cwh = [0, 2*vsize, vsize] + ind;  	% resorts so that first Nph cols are cens, then wids, then heights
		jacout(:, ind_cwh) = splitjac;

	else			% no jacobian
		splitgauss = calc_splitgauss_JAC_fast(eng,cen(ind),ht(ind),w1(ind),w2(ind));
		modelout = modelout + splitgauss';
	end
end


% tack on the backgrounds to the convolved splitgausses
modelout = modelout + constant;		% constant BG
modelout = modelout + (slope*eng);	% linear BG

