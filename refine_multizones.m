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


debug = 1;

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
%disp(' ');
%disp(' ');



%% === fit data ===

  opts = optimset ( ...
		'Jacobian',			'on', ...
		'MaxIter',			500,...
		'Display', 			'iter', ...
		'TolFun',			1e-16,
		'TolX',				1e-8,
		 );

	
	%% Check that array sizes 
	[~,maxsize] = computer;
	if system_octave;
		% https://github.com/calaba/octave-3.8.2-enable-64-ubuntu-14.04
		if sizemax < 2^31;
			warning(' Octave should probably be recompiled to handle larger arrays, see https://github.com/calaba/octave-3.8.2-enable-64-ubuntu-14.04');
		end
	end
	assert(length(varsin)*length(y_obs) < maxsize, 'The number of elements in the Jacobian is greater than the maxiumum array length that can be indexed.');


	if debug
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


if 1
	tic;
	disp('Beginning refinement...')
	[varsout,resnorm,resid,exitflag] = lsqnonlin( ...
										@(vars) objective(SYMS,vars,y_obs), ...
										varsin, bounds_lo, bounds_hi, opts);


	%% === post-process ===

		report_exitflag(exitflag);

		%	eliminates centers that had no data to fit
	%	varsout(find(SYMS{1}.VARS.freevars([1:end-1],1,1) == 0)) = NaN;
	disp(['Refinement complete!  Took ' num2str(i_obj) ' iterations.'])
	toc;

else
	varsout = varsin;
end

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
%	disp('Calling objective...')
	global i_obj
	if mod(i_obj, 10) == 0
		disp(['Objective iteration : ' num2str(i_obj)]);
		if mod(i_obj, 10) == 0
			toc;
		end
		if system_octave; fflush(stdout); end;
	end
	i_obj = i_obj + 1;

	weights = 0;		% DerivativeCheck passed when using weights
	func_mask = SYMS{1}.VARS.func_mask;

	if nargout > 1
		[func_out,jac_out] = calc_full_model(SYMS, vars);
		J = jac_out(func_mask, SYMS{1}.VARS.indfree);
	else
		func_out = calc_full_model(SYMS, vars);
	end

	if weights
		F = (func_out - y_obs) .* sqrt(SYMS{1}.VARS.w_obs);
	else
		F = func_out - y_obs;
	end
	F = F(func_mask);
end
end
