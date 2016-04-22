function SYM=generate_AUX(SYM,startvars);
% SYM=generate_AUX(SYM,startvars);
%	Each subset of data (a single DAT) has an AUX structure associated with it.
%	AUX contains auxiliary information about that DAT.
%
%	startvars(Nph x 2) has [cen,wid] for each mode
% 	Overall structure of AUX.auxvars:
%		Page 1:	[cens(Nph)		heights(Nph x Nq)...
%					0			bgConst			];
%
%		Page 2:	[widPh(Nph)		widRes(Nph x Nq)...
%					0			bgLinear		];
%


for ind=1:length(SYM)
	clear AUX;
	DAT=SYM{ind}.DAT;

	AUX.Nq=size(DAT.ydat,2);
	AUX.Nph=size(startvars,1);
	AUX.wdat=1./DAT.edat.^2;

	%new mask function makes mask and starting values for height fitting/decides whether to fit a height
	[AUX.mask goodheight free_cenht] = make_mask(SYM,startvars,ind);

	AUX.indE=[0 cumsum(sum(AUX.mask,1))];	% index of length of good points at each Q
	AUX.goodQ=sum(AUX.mask);


	goodcen=startvars(:,1);
	goodwid=zeros(size(goodcen));


	constantBackground = repmat(0,1,AUX.Nq);%min(DAT.ydat);

	reswids=merchop(SYM{ind}.Ei, SYM{ind}.chopfreq, goodcen);
	reswids=repmat(reswids(:),1,AUX.Nq);

	linearBackground = repmat(0,1,AUX.Nq);

	%%% AUX.auxvars
	%%% page1 is cen/height, page2 is phWid/resWid (resWid is assumed fixed, phWid could be fit)
	AUX.auxvars(:,:,1)=[goodcen goodheight;...
			    			0 constantBackground];
	AUX.auxvars(:,:,2)=[goodwid reswids;...
			    			0 linearBackground];

	% set BOUNDS
	AUX.bounds_L=zeros(size(AUX.auxvars));
	AUX.bounds_L([1:end-1],1,1)=min(DAT.xdat(:,1));		% min center

	AUX.bounds_H=ones(size(AUX.auxvars));
	AUX.bounds_H([1:end-1],1,1)=max(DAT.xdat(:,1));		% max cen is max of range
	AUX.bounds_H(:,[2:end],1)=1.5*max(max(DAT.ydat));	% maxheight is 1.5x maxdata
	AUX.bounds_H(:,:,2)=Inf;								% no max wid
	AUX.bounds_H(end,[2:end],2) = 1.5*max(max(DAT.ydat));

	% setup freevars field
	AUX.freevars=ones(size(AUX.auxvars));
	AUX.freevars(end,1,:) = 0;						% dummy variables
	AUX.freevars([1:end-1],:,1) = free_cenht;		% which heights are fitted is determined in make_mask.m
	AUX.freevars([1:end-1],[2:end],2)=0;			% res widths aren't free
	AUX.freevars([1:end-1],1,2)=1;					% phonon widths (0/1 fixed/free)
	AUX.freevars(end,[2:end],1) = 1;				% constant background (0/1 fixed/free)
	AUX.freevars(end,[2:end],2) = 1;				% linear background (0/1 fixed/free)
	AUX.indfree=find(AUX.freevars);

	AUX.peak_asymmetry=1.7;	%determined empirically for ARCS, using fityk

	% attach
	SYM{ind}.AUX=AUX;
end
