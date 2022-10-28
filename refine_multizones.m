function SYMS = refine_multizones(SYMS);

%% had to change to MATLAB syntax in these files:
#	modified:   0_subs/calc_full_model.m
#	modified:   0_subs/calc_singleQ.m
#	modified:   0_subs/make_aux_mask.m
#	modified:   assemble_VARS.m
#	modified:   make_full_jacobian.m
#	modified:   make_vars_mask.m
#	modified:   simulate_phonon_predictions.m

if system_octave;
    pkg load optim
end

global i_obj		% for counting how many iterations the objective function has executed
i_obj = 0;


DEBUG = 0;

VARS = SYMS{1}.VARS;
varsin = VARS.varsin(:);		% NOTE: varsin is only refineable variables

y_obs = VARS.y_obs;
w_obs = VARS.w_obs;

%display number of Q-points used
nZones = 0;
for ind = 1:length(SYMS)
	nZones = nZones + length(find(sum(SYMS{ind}.AUX.mask)));
end
disp(['  Refining with ' num2str(nZones) ' seperate Q-points']);

if isfield(SYMS{1}, 'CONFIG')
	CONFIG = SYMS{1}.CONFIG;
else
	CONFIG = {};
end


%% set optimization parameters
if isfield(CONFIG, 'optimset')
	disp(['Using optimization parameters defined in config.'])
else
	disp(['Using default optimization parameters'])
	opts = optimset ( ...
			'Jacobian',			'on', ...
			'MaxIter',			500,...
			'Display', 			'iter', ...
			'TolFun',			1e-18, ...
			'TolX',				1e-12,
			'CheckGradients',	1
			);
end

	

	%% === fit data ===
	%% Check that array sizes 
	[~,maxsize] = computer;
	if system_octave;
		% https://github.com/calaba/octave-3.8.2-enable-64-ubuntu-14.04
		if sizemax < 2^31;
			warning(' Octave should probably be recompiled to handle larger arrays, see https://github.com/calaba/octave-3.8.2-enable-64-ubuntu-14.04');
		end
	end
	assert(length(varsin)*length(y_obs) < maxsize, 'The number of elements in the Jacobian is greater than the maxiumum array length that can be indexed.');


	if DEBUG
		disp([' size(vars):' num2str(size(varsin))]);
		disp([' size(ydat):' num2str(size(y_obs))]);
		disp('  OPTS : ')
		opts
	end

if 0
	bounds_lo = VARS.bndsLO(VARS.indfree);
	bounds_hi = VARS.bndsHI(VARS.indfree);
else
	bounds_lo = VARS.bounds_L(VARS.indfree);
	bounds_hi = VARS.bounds_H(VARS.indfree);
end

tic;
disp('Beginning refinement...')
if 1
	disp("Using 'lsqnonlin'")
%	[varsout,resnorm,resid,exitflag] = lsqnonlin( ...
%										@(vars) objective(SYMS,vars,y_obs), ...
%										varsin, bounds_lo, bounds_hi, opts);

	[varsout,resnorm,resid,exitflag] = lsqnonlin( ...
										@(vars) objective(SYMS,vars,y_obs), ...
										varsin, bounds_lo, bounds_hi, opts);


	%% === post-process ===
	report_exitflag(exitflag);

elseif 1

	disp("Using 'lsqcurvefit'")
	[varsout,resnorm,resid,exitflag] = lsqcurvefit( ...
										@(vars) objective_curvefit(SYMS,vars,y_obs), ...
										varsin, bounds_lo, bounds_hi, opts);


	%% === post-process ===
	report_exitflag(exitflag);

else
	disp("Skipping refinement")
	varsout = varsin;
end
disp(['Refinement complete!  Took ' num2str(i_obj) ' iterations.'])
toc;

%% === update and uncertainty ===
SYMS=update_AUX(SYMS,varsout);
SYMS=assemble_VARS(SYMS);

%% check for NaN
%nansum_multizone = sum(isnan(SYMS{1}.AUX.auxvars(:,1,1)))


if 0
	disp('Calculating parameter uncertainty...')
	SYMS = calc_unc(SYMS);
end


function [F,J,varargout] = objective(SYMS, vars, y_obs);

	use_jac = (nargout > 1);
	USE_WEIGHTS = 0;

	%% track current iteration
	global i_obj
	if mod(i_obj, 10) == 0
		disp(['Objective iteration : ' num2str(i_obj)]);
		if mod(i_obj, 10) == 0
			toc;
		end
		if system_octave; fflush(stdout); end;
	end
	i_obj = i_obj + 1;


	%% calculate objective function (and derivative)
	func_mask = SYMS{1}.VARS.func_mask;

	if use_jac
		[func_out, jac_out] = calc_full_model(SYMS, vars);
		J = jac_out(func_mask, SYMS{1}.VARS.indfree);
	else
		func_out = calc_full_model(SYMS, vars);
	end

	% apply weights
	if USE_WEIGHTS
		weights = sqrt(SYMS{1}.VARS.w_obs);
		F = (func_out - y_obs) .* weights;
		if use_jac
			J = J .* repmat(weights(func_mask), 1, size(J, 2));
		end
	else
		F = func_out - y_obs;
	end
	F = F(func_mask);

	if DEBUG
	disp([' func_out :', num2str(size(func_out))])
	disp(['func_mask :', num2str(size(func_mask))])
	disp([' y_obs : ' num2str(size(y_obs))])
	disp([" F : ", num2str(size(F))])
		size(jac_out)
		tmp  = sum(jac_out);
		tmp(1:7)
		pause()
	end
end


function [F,J,varargout] = objective_curvefit(SYMS, vars, y_obs);

	use_jac = (nargout > 1);
	USE_WEIGHTS = 1;

	%% track current iteration
	global i_obj
	if mod(i_obj, 10) == 0
		disp(['Objective iteration : ' num2str(i_obj)]);
		if mod(i_obj, 10) == 0
			toc;
		end
		if system_octave; fflush(stdout); end;
	end
	i_obj = i_obj + 1;


	%% calculate objective function (and derivative)
	func_mask = SYMS{1}.VARS.func_mask;

	if use_jac
		[func_out,jac_out] = calc_full_model(SYMS, vars);
		J = jac_out(func_mask, SYMS{1}.VARS.indfree);
	else
		func_out = calc_full_model(SYMS, vars);
	end

	% apply weights
	if USE_WEIGHTS
%		weights = sqrt(SYMS{1}.VARS.w_obs);
		weights = SYMS{1}.VARS.w_obs;
		F = func_out .* weights;
		if use_jac
			J = J .* repmat(weights(func_mask), 1, size(J, 2));
		end
	end

	%
	F = F(func_mask);

	if DEBUG
	disp([' func_out :', num2str(size(func_out))])
	disp(['func_mask :', num2str(size(func_mask))])
	disp([' y_obs : ' num2str(size(y_obs))])
	disp([" F : ", num2str(size(F))])
		size(jac_out)
		tmp  = sum(jac_out);
		tmp(1:7)
		pause()
	end
end


end
