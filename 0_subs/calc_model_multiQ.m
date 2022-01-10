function [funcout,jacout]=calc_model_multiQ(SYM,varsin,varargin)
% [funcout,jout]=calc_model_multiQ_JAC(SYM,varsin,varargin)
%	Calculates the function and Jacobian (ensemble of partial derivatives) for
%	all Q in the SYM set.
%
%	Size of model for a single Q is just nEng.  Size of Jacobian for a single Q
%	is nEng x (3*nPhonon (cens,wids,hts) + BGconst + BGslope)
%
%	jacaux is the jacobian for all Q-points of a given SYM instance.
%		Size is nEng * nQ rows, cols is number of elements in AUX.auxvars

%(3*nPhonon + 2)*nQ cols. 
%		Each Q has ((cen-ht-wid)*nPhonon, BG)
%		many of those values will be zero, so sparse is preferred.


debug = 0;

%% check for NaN
%nansum_input = sum(isnan(SYM{1}.AUX.auxvars(:,1,1)))

if nargin > 1
	SYM=update_AUX(SYM,varsin);
end
VARS=SYM{1}.VARS;
Nph=VARS.Nph;
funcout=[];

%% check for NaN
%nansum_aux = sum(isnan(SYM{1}.AUX.auxvars(:,1,1)))

% make sparse jacobian
jacout = sparse(length(VARS.ydatin(:)), length(VARS.varsin(:)));

for inds=1:length(SYM)
	ycalc = [];
	ind_wids=[];
	clear AUX;
	clear DAT;

	AUX=SYM{inds}.AUX;
	DAT=SYM{inds}.DAT;

	N_good=length(find(AUX.mask));
	jacaux=sparse(N_good, numel(AUX.auxvars));


	for ind=1:AUX.Nq

		[model,jacobian] = calc_JAC_singleQ(AUX, DAT, ind);

		%%% now begins the indexing %%%
		%% “Lasciate ogni speranza, voi ch’entrate!” %%

		% The objective function is a stack of nEng for each Q.  So for each Q,
		% this pulls the correct indices
		this_indE= [ AUX.indE(ind)+1 : AUX.indE(ind+1) ];		% lineup for eng
		Qmask = AUX.mask(:,ind);

		%% THIS IS FRAGILE.  Ideally would build mask from DAT.ydat.  Then energy
		%% indexing will be more transparent, and the correct rows could be pulled
		%% using find(mask).  Only issue is that requires generating the entire
		%% jacobian, which may run into memory limits.

		%% estimating jacobian size:
		%
		%	for gamma point:
		%		nEng = 220, nQ = 102, nPhonon = 60
		% 		jacsize = (220 * 102) x (102 * 60)
		%		jacsize ~ 132 E+6 elements (so 1 GB ram for double-precision)
		%
		%	for general point:
		%		nEng = 220, nQ = 5k (100 zones * 48), nPhonon = 240  (checked, actually 360 phonons, need to recompute)
		%		jacsize = (220x5k) x (5k * 240)
		%		jacsize ~ 1.3 E+12 elements (so 10 TB ram!!!!!)
		%
		%	Selecting only the active-point subset reduces the size of 
		%		nEng 		by 3 (7k valid energy points out of 102 Q)
		%		nPhonon 	by (1.5 possibly optimistic)

		%		Gamma-point:
		%		nEng = 70, nQ = 100, nPhonon = 40
		%		jacsize = (70*100) x (100*40)
		%		jacsize ~ 28 E+6

		%		General-symmetry (nQ is 100 zones * 48) : 
		%		jacsize = (70*100*48) x (100*48*200)
		%		jacsize ~ 322 E+9 elements (2.5 TB for doubles)

		%	Generally, jacobian size goes like nEng * nPhonon * (nSym * nQ)**2
		%	going to sparse matrix reduces this even more.  Most of the array is
		%	heights, so going to spares reduces size by a bit under a factor of
		%	(nSym * nQ), or roughly 5,000.  Brings size back to a few GB


		jacobian(~isfinite(jacobian))=0;%hopefully this isn't necessary
		%model(~isfinite(model))=0;%hopefully this isn't necessary
		ycalc = [ycalc model(Qmask)];
		if ~isempty(find(isnan(jacobian)))
			disp(['Jacobian NaN at ind=' num2str(ind) ', SYM ' num2str(inds)]);
		end
		if ~isempty(find(isinf(jacobian)))
			disp(['Jacobian Inf at ind=' num2str(ind) ', SYM ' num2str(inds)]);
		end
		% place subsect of jacobian into correct part of jacaux

		assert(isequal(length(this_indE),length(find(Qmask))));

		indcen=[1:Nph];
		indwid=[1:Nph]+((Nph+1)*(AUX.Nq+1));	% Nph wids, starting after hts/cens
		indht =[1:Nph]+ ind*(Nph+1);
		indconstBG =(ind+1)*(Nph+1);
		indlinBG =(Nph+1)*(2+AUX.Nq+ind);


		% chunks of heights, each is Nph cols
		assert(isequal(size(jacaux(this_indE,indcen)),size(jacobian(Qmask,3*[1:Nph]-2))));
		jacaux(this_indE,indcen) = jacobian(Qmask,3*[1:Nph]-2);	% centers
		jacaux(this_indE, indht) = jacobian(Qmask,3*[1:Nph]-1);	% heights

		% widths may be wrong - 1st col should be intrinsic, others should be resfunc (so fixed, and never included)
		jacaux(this_indE,indwid) =jacobian(Qmask,3*[1:Nph]);	% widths (starting halfway)	-why is factor of 0.5 needed to agree with finite-diff?

		% background terms for the jacobian
		jacaux(this_indE,indconstBG) =jacobian(Qmask,end-1);
		jacaux(this_indE,indlinBG) =jacobian(Qmask,end);
	end


	ind_wids=[ind_wids indwid(1)];
	SYM{inds}.funaux=ycalc';
	SYM{inds}.jacaux=jacaux;
	funcout=[funcout; ycalc'];%funcout=[funcout; ycalc(AUX.mask)];

	if 0
		disp(['indfree: ' num2str(length(AUX.indfree))]);
		disp([' ycalc : ' num2str(length(AUX.mask))]);
		disp([' jacaux : ' num2str(size(jacaux))]);
		disp([' jacout : ' num2str(size(jout))]);
	end
end

funcout=funcout(:);
% now combine jacobians
%  This could be cleaned up a bit to improve the speed, see:
%	http://www.mathworks.com/help/matlab/math/accessing-sparse-matrices.html
%
%	Basically, it's better to build index arrays for row, col, and then run
%	sparse( iRows, iCols, values)
%
%	The following section *almost* does that


ind_func=0;
ind_qs=0;
%ind_wids;
Nqs=[0 VARS.Nqs];
for ind=1:length(SYM)
	%row of jacout
	ind_func=[ind_func length(find(SYM{ind}.AUX.mask))];
	jrow=[ind_func(ind)+1 : sum(ind_func)];

	%columns of jacaux
	indcen = [1:Nph];
	indwid = [1:Nph]+((Nph+1)*(VARS.Nqs(ind)+1));
	indht = [];
	indconst = [];
	indlin = [];
	for idx = 1:VARS.Nqs(ind)
		indht = [indht [1:Nph]+idx*(Nph+1)];
		indconst = [indconst (idx+1)*(Nph+1)];
		indlin = [indlin (Nph+1)*(VARS.Nqs(ind)+2+idx)];
	end

	%columns of jacout
	idxcen = [1:Nph];
	idxwid = [1:Nph]+(Nph+1)*(1+sum(VARS.Nqs));
	idxht = indht + (Nph+1)*sum(Nqs([1:ind]));
	idxconst = indconst + (Nph+1)*sum(Nqs([1:ind]));
	idxlin = idxconst + (Nph+1)*(1+sum(VARS.Nqs));

	%assigns values

	if debug
		disp(' ')
		ind
		size(jacout)
		disp([' jacaux : ' num2str(size(SYM{ind}.jacaux)) ])
		length(jrow)
		idxcen
		indcen
	end

	jacout(jrow,idxcen)=SYM{ind}.jacaux(:,indcen);			% centers
	jacout(jrow,idxwid)=SYM{ind}.jacaux(:,indwid);			% widths
	jacout(jrow,idxht)=SYM{ind}.jacaux(:,indht);			% heights
	jacout(jrow,idxconst)=SYM{ind}.jacaux(:,indconst);		% constant BG
	jacout(jrow,idxlin)=SYM{ind}.jacaux(:,indlin);			% linear BG
end
jacout=jacout(:,SYM{1}.VARS.indfree);



%display each iteration of fitting
if 0;
	eng = SYM{1}.DAT.xdat(AUX.mask);
	%semilogy(SYM{1}.DAT.xdat(AUX.mask),SYM{1}.DAT.ydat(AUX.mask),'bo',SYM{1}.DAT.xdat(AUX.mask),funcout,'ro')
	hold off; errorbar(eng, SYM{1}.DAT.ydat(AUX.mask), SYM{1}.DAT.edat(AUX.mask),'b--');
	length(funcout)
	length(eng)
	hold on; plot(eng, funcout, 'r-', 'linewidth',1);
	pause(.01)
end
%toc

if system_octave; fflush(stdout); end
