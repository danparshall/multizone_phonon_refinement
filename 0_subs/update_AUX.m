function SYMS=update_AUX(SYMS,varsin);


VARS=SYMS{1}.VARS;
VARS.allvars(VARS.indfree)=varsin;



for ind=1:length(SYMS)
	if 0  % old indexing
	% update AUX.auxvars
	Nqs=[0 VARS.Nqs];	% leading 0 for indexing  VARS.Nqs is is number in each AUXVARS
	cNqs=cumsum(Nqs);
		SYMS{ind}.AUX.auxvars(:,1,:) = VARS.allvars(:,1,:);
		SYMS{ind}.AUX.auxvars(:, 1+[1:Nqs(ind+1)], :) = VARS.allvars(:, 1+[cNqs(ind)+1:cNqs(ind+1)], :);
		if length([1:Nqs(ind+1)]) ~= length([cNqs(ind)+1:cNqs(ind+1)]); disp(' oops'); end
	else
		freevars = SYMS{ind}.AUX.freevars;
%		auxvars = SYMS{ind}.AUX.auxvars;
		vars_mask = SYMS{ind}.AUX.vars_mask;
		SYMS{ind}.AUX.auxvars(find(freevars)) = VARS.allvars(vars_mask(find(freevars)));
	end
end
SYMS{1}.VARS=VARS;
