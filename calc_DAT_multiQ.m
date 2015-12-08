function [funcout,jacout]=calc_DAT_multiQ(SYMDAT,varsin)
% [funcout,jout]=calc_model_multiQ_JAC(SYMDAT,varsin)
%	Calculates the function and Jacobian (ensemble of partial derivatives) for
%	all Q in the SYMDAT set.

SYMDAT=update_AUX(SYMDAT,varsin);
VARS=SYMDAT{1}.VARS;
Nph=VARS.Nph;
funcout=[];
jacout=[];

for inds=1:length(SYMDAT)
	ycalc = [];
	ind_wids=[];
	clear AUX;
	clear SYM;

	AUX=SYMDAT{inds}.AUX;
	SYM=SYMDAT{inds}.SQW;

	%ycalc=zeros(size(SYM.ydat));
	N_good=length(find(AUX.mask));
	%ycalc=zeros(size(SYM.xdat,1),AUX.Nq);
	jacux=zeros(N_good, numel(AUX.auxvars));
%	jacaux=sparse( [],[],[], N_good, numel(AUX.auxvars));
%	jacaux=full(jacaux);

	for ind=1:AUX.Nq

	% energy used in singleQ comes from mask(:,ind)
		%if confirm_dataset_valid(SYM.xdat(AUX.mask(:,ind),ind));
			[modout,thisjac]=calc_JAC_singleQ(AUX,SYM,ind);
		%else
		%	modout = zeros([1 0]);
		%	thisjac = zeros([0 3]);
		%end

		thisjac(~isfinite(thisjac))=0;%hopefully this isn't necessary
		%modout(~isfinite(modout))=0;%hopefully this isn't necessary
		ycalc = [ycalc modout];
		if ~isempty(find(isnan(thisjac)))
			disp(['Jacobian NaN at ind=' num2str(ind) ', SYMDAT ' num2str(inds)]);
		end
		if ~isempty(find(isinf(thisjac)))
			disp(['Jacobian Inf at ind=' num2str(ind) ', SYMDAT ' num2str(inds)]);
		end
		% place subsect of thisjac into correct part of jacaux
		this_indE= [ AUX.indE(ind)+1 : AUX.indE(ind+1) ];		% lineup for eng
		if length(this_indE)~=length(find(AUX.mask(:,ind)));	disp(' ERROR in "calc_DAT_multi_Q"'); end

		indcen=[1:Nph];
		indwid=[1:Nph]+((Nph+1)*(AUX.Nq+1));	% Nph wids, starting after hts/cens
		indht =[1:Nph]+ ind*(Nph+1);
		indconstBG =(ind+1)*(Nph+1);
		indlinBG =(Nph+1)*(2+AUX.Nq+ind);


		% chunks of heights, each is Nph long
		jacaux(this_indE,indcen) =thisjac(:,3*[1:Nph]-2);	% centers
		jacaux(this_indE, indht) =thisjac(:,3*[1:Nph]-1);	% heights

		% widths may be wrong - 1st col should be intrinsic, others should be resfunc (so fixed, and never included)
		jacaux(this_indE,indwid) =thisjac(:,3*[1:Nph]);	% widths (starting halfway)	-why is factor of 0.5 needed to agree with finite-diff?

		% background terms for the jacobian
		jacaux(this_indE,indconstBG) =thisjac(:,end-1);
		jacaux(this_indE,indlinBG) =thisjac(:,end);
	end


	ind_wids=[ind_wids indwid(1)];
	SYMDAT{inds}.funaux=ycalc';
	SYMDAT{inds}.jacaux=jacaux;
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
jacout=zeros( length(funcout), numel(VARS.allvars));
%jacout=sparse( [],[],[], length(funcout), numel(VARS.allvars));
%jacout=full(jacout);
ind_func=0;
ind_qs=0;
%ind_wids;
Nqs=[0 VARS.Nqs];
for ind=1:length(SYMDAT)
	%row of jacout
	ind_func=[ind_func length(find(SYMDAT{ind}.AUX.mask))];
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
	jacout(jrow,idxcen)=SYMDAT{ind}.jacaux(:,indcen);			% centers
	jacout(jrow,idxwid)=SYMDAT{ind}.jacaux(:,indwid);			% widths
	jacout(jrow,idxht)=SYMDAT{ind}.jacaux(:,indht);			% heights
	jacout(jrow,idxconst)=SYMDAT{ind}.jacaux(:,indconst);		% constant BG
	jacout(jrow,idxlin)=SYMDAT{ind}.jacaux(:,indlin);			% linear BG
end
jacout=jacout(:,SYMDAT{1}.VARS.indfree);
%display each iteration of fitting
if 0;
	eng = SYMDAT{1}.SQW.xdat(AUX.mask);
	%semilogy(SYMDAT{1}.SQW.xdat(AUX.mask),SYMDAT{1}.SQW.ydat(AUX.mask),'bo',SYMDAT{1}.SQW.xdat(AUX.mask),funcout,'ro')
	hold off; errorbar(eng, SYMDAT{1}.SQW.ydat(AUX.mask), SYMDAT{1}.SQW.edat(AUX.mask),'b--');
	length(funcout)
	length(eng)
	hold on; plot(eng, funcout, 'r-', 'linewidth',1)
	pause(.01)
end
%toc
fflush(stdout);
