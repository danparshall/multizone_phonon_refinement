function SYMS=update_AUX(SYMS,varsin);
% update values from SYMS.allvars (accessed by leasqr) to various AUX.auxvars (used by "calc_DAT")

%SYMS{1}.VARS.allvars
% update VARS.allvars
VARS=SYMS{1}.VARS;
VARS.allvars(VARS.indfree)=varsin;

% update AUX.auxvars
Nqs=[0 VARS.Nqs];	% leading 0 for indexing
cNqs=cumsum(Nqs);

for ind=1:length(SYMS)
	SYMS{ind}.AUX.auxvars(:,1,:)=VARS.allvars(:,1,:);
%1+[1:Nqs(ind+1)]
%1+[cNqs(ind)+1:cNqs(ind+1)]
	SYMS{ind}.AUX.auxvars(:, 1+[1:Nqs(ind+1)], :)=VARS.allvars(:,1+[cNqs(ind)+1:cNqs(ind+1)],:);
	if length([1:Nqs(ind+1)])~=length([cNqs(ind)+1:cNqs(ind+1)]); disp(' oops');end

chk=SYMS{ind}.AUX.auxvars(:,:,1);
%pause
end
SYMS{1}.VARS=VARS;
