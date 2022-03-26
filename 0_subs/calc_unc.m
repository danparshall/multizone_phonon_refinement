function SYMS = calc_unc(SYMS)
% use QR decomposition to estimate uncertainty of parameters


ts = time;
disp(['  Estimating parameter uncertainty...'])
debug = 0;

VARS = SYMS{1}.VARS;
func_mask = VARS.func_mask;

y_obs = VARS.y_obs(func_mask);
w_obs = VARS.w_obs(func_mask);		% weights

if debug;
    tot_pts = [];
    tot_valid = [];
    for i_sym = 1:length(SYMS)
        tot_pts = [tot_pts, prod(size(SYMS{i_sym}.DAT.edat)) ];
        tot_valid = [tot_valid, length(find(SYMS{i_sym}.DAT.edat > 0)) ];
    end

    tot_pts
    tot_valid
    length(VARS.w_obs)

    length(w_obs)
    sum(isnan(w_obs))
end

[func, jac] = calc_full_model(SYMS, VARS.varsin);

residuals = func(func_mask) - y_obs;
%resnorm = residuals.^2 .* w_obs;
resnorm = residuals.^2;

jac = jac(func_mask, VARS.indfree);
[N_pts,N_vars] = size(jac);
disp('  .... decomposing jacobian')
[Q,R] = qr(jac,0);
covar = sum(inv(R).^2, 2);

% uncertainty of each parameter of varsin
VARS.uncertainty(VARS.indfree) = full(sqrt(covar*sum(resnorm)/(N_pts-N_vars)));

% update the uncertainty for each parameter in AUX.auxvars
for ind=1:length(SYMS)
	freevars = SYMS{ind}.AUX.freevars;
	vars_mask = SYMS{ind}.AUX.vars_mask;
    SYMS{ind}.AUX.uncertainty = zeros(size(freevars));
	SYMS{ind}.AUX.uncertainty(find(freevars)) = VARS.uncertainty(vars_mask(find(freevars)));
end

te = time;
tdelta = te - ts;
disp(['    .... finished; elapsed time ' num2str(tdelta)])