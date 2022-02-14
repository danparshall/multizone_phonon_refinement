function SYMS = update_AUX(SYMS,varsin);

VARS = SYMS{1}.VARS;
VARS.allvars(VARS.indfree) = varsin;

for ind=1:length(SYMS)
	freevars = SYMS{ind}.AUX.freevars;
	vars_mask = SYMS{ind}.AUX.vars_mask;
	% the positon in vars_mask corresponds to auxvars; the value in vars_mask is the index (column of jac) that variable occupies in VARS.allvars
	SYMS{ind}.AUX.auxvars(find(freevars)) = VARS.allvars(vars_mask(find(freevars)));
end

SYMS{1}.VARS = VARS;