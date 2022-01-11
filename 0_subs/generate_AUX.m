function SYMS=generate_AUX(SYMS);
% SYMS=generate_AUX(SYMS);
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
%  [13,1] = peak_aSYMSmetry


for ind = 1:length(SYMS)
	clear AUX;
	DAT = SYMS{ind}.DAT;
	startvars = SYMS{ind}.startvars;

	AUX.Nq = size(DAT.ydat,2);
	AUX.Nph = size(startvars,1);
	AUX.wdat = 1./DAT.edat.^2;

	if isfield(DAT,'eng');
		AUX.eng = DAT.eng;
	else;
		AUX.eng = DAT.xdat(:,1);
	end

	% makes mask and starting values for height fitting, and decides whether to fit a height
	[AUX.mask goodheight free_cenht] = make_mask(SYMS,startvars,ind);
	free_cenht
	AUX.indE = [0 cumsum(sum(AUX.mask,1))];	% index of length of good points at each Q
	AUX.goodQ = sum(AUX.mask);


	goodcen = startvars(:,1);
	goodwid = startvars(:,2);
	goodheight = startvars(:,3:end);


	constantBackground = repmat(0,1,AUX.Nq); %min(DAT.ydat);

	% NOTE : we are define the resolution width when we initialize.  In principle we should update the resolution width,
	% based on our best guess about the peak center; but getting that into the Jacobian would be really difficult.
	% As a compromise, we just leave the reswidth fixed (not too big a deal, if our original DFT cens weren't far off).
	reswids = merchop(SYMS{ind}.Ei, SYMS{ind}.chopfreq, goodcen);
	reswids = repmat(reswids(:),1,AUX.Nq);

	linearBackground = repmat(0,1,AUX.Nq);

	%%% AUX.auxvars
	%%% page1 is cen/height, page2 is phWid/resWid (resWid is assumed fixed, phWid could be fit)
	AUX.auxvars(:,:,1) = [goodcen goodheight;...
			    			0 constantBackground];
	AUX.auxvars(:,:,2) = [goodwid reswids;...
			    			0 linearBackground];

	% set BOUNDS (have the same structure as auxvars)
	AUX.bounds_L = 0.01 * ones(size(AUX.auxvars));
	AUX.bounds_L([1:end-1],1,1) = min(DAT.xdat(:,1));		% min center

	AUX.bounds_H = ones(size(AUX.auxvars));
%	AUX.bounds_H([1:end-1],1,1) = max(DAT.xdat(:,1));		% max cen is max of range
	AUX.bounds_H([1:end-1],1,1) = 10*goodcen;				% max cen is 10x the DFT prediction
	AUX.bounds_H(:,[2:end],1) = 1.5*max(max(DAT.ydat));		% maxheight is 1.5x maxdata
%	AUX.bounds_H(:,:,2) = Inf;								% no max wid

	AUX.bounds_H([1:end-1],[2:end],2) = repmat(10*reswids(:,1), 1, AUX.Nq);					% wids, set to 10x the resolution
	AUX.bounds_H(end,[2:end],1) = max(max(DAT.ydat));		% constant BG max set to max of observed counts

	delta_counts_each_Q = max(DAT.ydat) - min(DAT.ydat);
	delta_energy = max(AUX.eng) - min(AUX.eng);
	AUX.bounds_H(end,[2:end],2) = delta_counts_each_Q/delta_energy;	% linear BG max

	% setup freevars field (same structure as auxvars, defines which params can be fit)
	AUX.freevars = ones(size(AUX.auxvars));
	AUX.freevars(end,1,:) = 0;						% placeholder variables, not actually fit
	AUX.freevars([1:end-1],:,1) = free_cenht;		% which heights are fitted is determined in make_mask.m
	AUX.freevars([1:end-1],[2:end],2)=0;			% res widths aren't free
	AUX.freevars([1:end-1],1,2) = 1;				% phonon widths (0/1 fixed/free)
	AUX.freevars(end,[2:end],1) = 1;				% constant background (0/1 fixed/free)
	AUX.freevars(end,[2:end],2) = 1;				% linear background (0/1 fixed/free)
	AUX.indfree = find(AUX.freevars);

	AUX.peak_asymmetry = 1.7;		% How wide is low-energy side, relative to high-energy side.  Determined empirically for ARCS, using fityk.

	% attach
	SYMS{ind}.AUX = AUX;
end
