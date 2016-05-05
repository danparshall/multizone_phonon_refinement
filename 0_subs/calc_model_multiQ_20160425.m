function [funcout,jacout]=calc_DAT_multiQ(SYM,varsin)
% [funcout,jout]=calc_model_multiQ_JAC(SYM,varsin)
%	Calculates the function and Jacobian (ensemble of partial derivatives) for
%	all Q in the SYM set.

debug = 0;

SYM=update_AUX(SYM,varsin);
VARS=SYM{1}.VARS;
Nph=VARS.Nph;
funcout=[];
jacout=[];

for inds=1:length(SYM)
	ycalc = [];
	ind_wids=[];
	clear AUX;
	clear DAT;

	AUX=SYM{inds}.AUX;
	DAT=SYM{inds}.DAT;

	%ycalc=zeros(size(DAT.ydat));
	N_good=length(find(AUX.mask));
	%ycalc=zeros(size(DAT.xdat,1),AUX.Nq);
	jacaux=zeros(N_good, numel(AUX.auxvars));
%	jacaux=sparse( [],[],[], N_good, numel(AUX.auxvars));
%	jacaux=full(jacaux);

	for ind=1:AUX.Nq

	% energy used in singleQ comes from mask(:,ind)
		%if confirm_dataset_valid(DAT.xdat(AUX.mask(:,ind),ind));
			[model,jacobian]=calc_JAC_singleQ(AUX,DAT,ind);
		%else
		%	model = zeros([1 0]);
		%	jacobian = zeros([0 3]);
		%end

		jacobian(~isfinite(jacobian))=0;%hopefully this isn't necessary
		%model(~isfinite(model))=0;%hopefully this isn't necessary
		ycalc = [ycalc model];
		if ~isempty(find(isnan(jacobian)))
			disp(['Jacobian NaN at ind=' num2str(ind) ', SYM ' num2str(inds)]);
		end
		if ~isempty(find(isinf(jacobian)))
			disp(['Jacobian Inf at ind=' num2str(ind) ', SYM ' num2str(inds)]);
		end
		% place subsect of jacobian into correct part of jacaux
		this_indE= [ AUX.indE(ind)+1 : AUX.indE(ind+1) ];		% lineup for eng

		assert(length(this_indE)==length(find(AUX.mask(:,ind))));

		indcen=[1:Nph];
		indwid=[1:Nph]+((Nph+1)*(AUX.Nq+1));	% Nph wids, starting after hts/cens
		indht =[1:Nph]+ ind*(Nph+1);
		indconstBG =(ind+1)*(Nph+1);
		indlinBG =(Nph+1)*(2+AUX.Nq+ind);


		% chunks of heights, each is Nph long
size(jacaux(this_indE,indcen))
size(jacobian(:,3*[1:Nph]-2))
		assert(size(jacaux(this_indE,indcen)) == size(jacobian(:,3*[1:Nph]-2)));
		jacaux(this_indE,indcen) =jacobian(:,3*[1:Nph]-2);	% centers
		jacaux(this_indE, indht) =jacobian(:,3*[1:Nph]-1);	% heights

		% widths may be wrong - 1st col should be intrinsic, others should be resfunc (so fixed, and never included)
		jacaux(this_indE,indwid) =jacobian(:,3*[1:Nph]);	% widths (starting halfway)	-why is factor of 0.5 needed to agree with finite-diff?

		% background terms for the jacobian
		jacaux(this_indE,indconstBG) =jacobian(:,end-1);
		jacaux(this_indE,indlinBG) =jacobian(:,end);
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
jacout=zeros( length(funcout), numel(VARS.allvars));
%jacout=sparse( [],[],[], length(funcout), numel(VARS.allvars));
%jacout=full(jacout);
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
	hold on; plot(eng, funcout, 'r-', 'linewidth',1)
	pause(.01)
end
%toc
fflush(stdout);
