function SYM=update_AUX(SYM,varsin);
% update values from VARS.allvars (accessed by lsqnonlin) to various AUX.auxvars (used by "calc_model")


% update VARS.allvars
VARS=SYM{1}.VARS;
VARS.allvars(VARS.indfree)=varsin;

% update AUX.auxvars
Nqs=[0 VARS.Nqs];	% leading 0 for indexing  VARS.Nqs is is number in each AUXVARS

cNqs=cumsum(Nqs);

for ind=1:length(SYM)

	% update centers
	SYM{ind}.AUX.auxvars(:,1,:)=VARS.allvars(:,1,:);


	% update everything else
	assert(length([1:Nqs(ind+1)])==length([cNqs(ind)+1:cNqs(ind+1)]), ' index error in update_AUX');
	SYM{ind}.AUX.auxvars(:, 1+[1:Nqs(ind+1)], :)=VARS.allvars(:,1+[cNqs(ind)+1:cNqs(ind+1)],:);

end

SYM{1}.VARS=VARS;
