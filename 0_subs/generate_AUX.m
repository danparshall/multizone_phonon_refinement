function SYMS = generate_AUX(SYMS);
% SYMS = generate_AUX(SYMS);
%	Each subset of data (a single DAT) has an AUX structure associated with it.
%	AUX contains auxiliary information about that DAT.
%
%	startvars is [cens wids heights], size = nPhonon x (nQ+2)
%
%
% 	Overall structure of AUX.auxvars:
%		Page 1:	[cens(Nph)		heights(Nph x Nq)...
%					0			bgConst			];
%
%		Page 2:	[widPh(Nph)		widRes(Nph x Nq)...
%					0			bgLinear		];
%
%
%	Size estimate for BKBO off-symm:
%		nCen = [6 * 60]
%		nQ = 48 * 100
%		size ~ 2 * nCen * nQ * 8 bytes = 28 MB

%% fields are:
%  [1,1] = Nq
%  [2,1] = Nph
%  [3,1] = wdat
%  [4,1] = eng
%  [5,1] = mask
%  [6,1] = indE
%  [7,1] = goodQ
%  [8,1] = auxvars
%  [9,1] = bounds_L
%  [10,1] = bounds_H
%  [11,1] = freevars
%  [12,1] = indfree
%  [13,1] = peak_asymmetry

debug = 1;
n_cens = size(SYMS{1}.startvars, 1);

for i_sym = 1:length(SYMS)
	disp('')
	disp(["STARTING AUX FOR SYM : " num2str(i_sym)])
	clear AUX;
	DAT = SYMS{i_sym}.DAT;
	startvars = SYMS{i_sym}.startvars;

	if isfield(DAT,'eng');
		AUX.eng = DAT.eng;
	else;
		AUX.eng = DAT.x_dat(:,1);
	end

	% makes mask and starting values for height fitting, and decides whether to fit a height
	% mask is same dimensions as data
	[AUX.mask goodheight free_cenht] = make_aux_mask(SYMS,startvars,i_sym);

	goodcen = startvars(:,1);
	goodwid = startvars(:,2);
	goodheight = startvars(:,3:end);


	% only allow this Q to be used if it has enough data points
	% Specifically, we require 2 for BG, plus 3 for every peak whose center is within range
	qpoints = sum(AUX.mask, 1);
	fittable_Q = find(qpoints > (2 + 3*sum(free_cenht(:, 2:end), 1) ));
	goodheight = goodheight(:, fittable_Q);
	free_cenht = free_cenht(:, [1, 1+fittable_Q]);
	AUX.mask = AUX.mask(:, fittable_Q);
	DAT.y_dat = DAT.y_dat(:, fittable_Q);
	DAT.e_dat = DAT.e_dat(:, fittable_Q);
	DAT.Q_hkl = DAT.Q_hkl(fittable_Q, :);
	DAT.Qmags = DAT.Qmags(fittable_Q);

	% now set Nq
	AUX.Nq = size(DAT.y_dat,2);
	AUX.Nph = size(startvars,1);

	% user-provided fit settings can be used, otherwise just take the defaults
	if isfield(SYMS{1}, 'CONFIG')
		CONFIG = SYMS{1}.CONFIG;
		if length(CONFIG.refine_centers) > 0
			assert(length(CONFIG.refine_centers) == AUX.Nph, "The number of phonons specified in CONFIG.refine_centers inconsistent with AUX.Nph")
			refine_cens = CONFIG.refine_centers;
			refine_wids = CONFIG.refine_widths;
		else
			refine_cens = free_cenht(:, 1);
			refine_wids = refine_cens;
		end
		refine_const_bg = CONFIG.refine_const_bg;
		refine_linear_bg = CONFIG.refine_linear_bg;
	else
		disp('No CONFIG file specified; using defaults for fitting centers, widths, and backgrounds.')
		refine_cens = free_cenht(:, 1);
		refine_wids = refine_cens;
		refine_const_bg = 1;
		refine_linear_bg = 0;
	end

	% NOTE : we define the resolution width when we initialize.  In principle we should update the resolution width
	% based on our current best guess about the peak center; but getting that into the Jacobian would be really difficult.
	% As a compromise, we just leave the reswidth fixed (not too big a deal, if our original DFT cens weren't far off).
	reswids = merchop(SYMS{i_sym}.DAT.Ei, SYMS{i_sym}.DAT.chopfreq, goodcen);
%	reswids = merchop(SYMS{i_sym}.Ei, SYMS{i_sym}.chopfreq, goodcen);
	reswids = repmat(reswids(:), 1, AUX.Nq);

	tmp = DAT.y_dat;
	tmp(~AUX.mask) = nan;
%	const_BG = min(tmp);
	const_BG = zeros(1,AUX.Nq);
%	const_BG = min(DAT.y_dat ./ AUX.mask)   % division by mask is to force the out-of-bounds areas to NaN, rather than zeros
	linear_BG = zeros(1,AUX.Nq);

	%%% AUX.auxvars
	%%% page1 is cen/height, page2 is phWid/resWid (resWid is assumed fixed, phWid can be fit)
	AUX.auxvars(:,:,1) = [goodcen	goodheight;...
			    			0 		const_BG];
	AUX.auxvars(:,:,2) = [goodwid 	reswids;...
			    			0 		linear_BG];

	% set BOUNDS (have the same structure as auxvars); 
	% NOTE: https://www.mathworks.com/help/optim/ug/iterations-can-violate-constraints.html
	delta_counts = max(max(DAT.y_dat(AUX.mask))) - min(min(DAT.y_dat(AUX.mask)))
	delta_energy = max(AUX.eng) - min(AUX.eng);

	AUX.bounds_L = 0.001 * ones(size(AUX.auxvars));
	AUX.bounds_L([1:end-1], 1, 1) = 0.1*min(goodcen);				% min center is 10% of minimum DFT prediction
	AUX.bounds_L([end, 2:end, 1]) = -0.1 * delta_counts;			% min constant BG set to -10% of observed intensity range
	AUX.bounds_L([end, 2:end, 2]) = -0.1*delta_counts/delta_energy;	% min linear BG


	AUX.bounds_H = ones(size(AUX.auxvars));
	AUX.bounds_H([1:end-1],1,1) = 10*goodcen;					% max cen is 10x the DFT prediction
	AUX.bounds_H(:,[2:end],1) = 1.5*max(max(DAT.y_dat));			% maxheight is 1.5x maxdata
	AUX.bounds_H(1:end-1,1,2) = 10*max(goodwid);				% max width is 10x the max DFT prediction
	AUX.bounds_H(end,[2:end],1) = 0.5*max(max(DAT.y_dat));		% constant BG max set to (max of observed counts)/2
	AUX.bounds_H(end,[2:end],2) = delta_counts/delta_energy;	% linear BG max
	AUX.bounds_H([1:end-1],[2:end],2) = repmat(10*abs(reswids(:,1)), 1, AUX.Nq);	% set to 10x the resolution (abs bc non-kinematic become negative)

	% setup freevars field (same structure as auxvars, defines which params can be fit)
	AUX.freevars = ones(size(AUX.auxvars));
	AUX.freevars(end,1,:) = 0;						% placeholder variables, not actually fit
	AUX.freevars([1:end-1],:,1) = free_cenht;		% which heights are fitted is determined in make_aux_mask.m
	AUX.freevars([1:end-1],[2:end],2)=0;			% res widths aren't free
	AUX.freevars(end,[2:end],1) = refine_const_bg;	% constant background (0/1 fixed/free)
	AUX.freevars(end,[2:end],2) = refine_linear_bg;	% linear background (0/1 fixed/free)

	% which centers/widths are refined is set either in CONFIG, or in aux_make_mask.
	AUX.freevars([1:end-1],1,2) = refine_wids;		% phonon widths (0/1 fixed/free)
	AUX.freevars([1:end-1],1,1) = refine_cens;		% phonon centers (using )
	AUX.indfree = find(AUX.freevars);

	AUX.peak_asymmetry = DAT.peak_asymmetry;
%	AUX.peak_asymmetry = SYMS{i_sym}.peak_asymmetry;


	% attach
	SYMS{i_sym}.AUX = AUX;
	SYMS{i_sym}.DAT = DAT;
end
