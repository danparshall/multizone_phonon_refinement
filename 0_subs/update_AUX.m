function SYM=update_AUX(SYM,varsin);
% update values from SYM.allvars (accessed by leasqr) to various AUX.auxvars (used by "calc_DAT")

%SYM{1}.VARS.allvars
% update VARS.allvars
VARS=SYM{1}.VARS;
VARS.allvars(VARS.indfree)=varsin;

% update AUX.auxvars
Nqs=[0 VARS.Nqs];	% leading 0 for indexing
cNqs=cumsum(Nqs);

for ind=1:length(SYM)
	SYM{ind}.AUX.auxvars(:,1,:)=VARS.allvars(:,1,:);
%1+[1:Nqs(ind+1)]
%1+[cNqs(ind)+1:cNqs(ind+1)]
	SYM{ind}.AUX.auxvars(:, 1+[1:Nqs(ind+1)], :)=VARS.allvars(:,1+[cNqs(ind)+1:cNqs(ind+1)],:);
	if length([1:Nqs(ind+1)])~=length([cNqs(ind)+1:cNqs(ind+1)]); disp(' oops');end

chk=SYM{ind}.AUX.auxvars(:,:,1);
%pause
end
SYM{1}.VARS=VARS;
