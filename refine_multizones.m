function SYMS = refine_multizones(SYMS);


debug = 1;

VARS=SYMS{1}.VARS;
varsin=VARS.varsin(:);

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
		'DerivativeCheck',	'off',...
		'MaxIter',			400,...
		'Display', 			'iter');

	
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
	[varsout,resnorm,resid,exitflag] = lsqnonlin( ...
										@(vars) objective(SYMS,vars,y_obs), ...
										varsin, bounds_lo, bounds_hi, opts);


	%% === post-process ===

		report_exitflag(exitflag);

		%	eliminates centers that had no data to fit
	%	varsout(find(SYMS{1}.VARS.freevars([1:end-1],1,1) == 0)) = NaN;
else
	varsout = varsin;
end

%% === update and uncertainty ===
SYMS=update_AUX(SYMS,varsout);
SYMS=assemble_VARS(SYMS);

%% check for NaN
%nansum_multizone = sum(isnan(SYMS{1}.AUX.auxvars(:,1,1)))


if 0
	unc = calc_unc(SYMS);
end



% === plotting ===
if 1
	funcout = calc_model_multiQ(SYMS,varsout);
	disp('Plotting...')
	for ind_sym = 1:length(SYMS)
		SYM = SYMS{ind_sym};

		disp(['  ind_sym : ' num2str(ind_sym)])
		disp('startvars :')
		SYM.startvars
		disp('Cens :')
		cen = SYMS{1}.VARS.allvars(1:end-1, 1, 1)
		funcindex = cumsum(sum(SYM.AUX.mask));
		mask = SYM.AUX.mask;
		DAT = SYM.DAT;
		Nq = SYM.AUX.Nq
		
		for ind = 1:Nq
			if sum(mask(:,ind))>0
				qpoint = DAT.HKL_vals(ind,:);
				xdata = DAT.xdat(mask(:,ind),ind);
				ydata = DAT.ydat(mask(:,ind),ind);
				edata = DAT.edat(mask(:,ind),ind);
%save_data(xdata',ydata',edata',qpoint);
%SYMS{1}.VARS.allvars(:,ind+1,:)
%SYMS{1}.AUX.freevars(:,ind+1,:)
				xfit = DAT.xdat(mask(:,ind),ind);
				if ind ~= 1
%					yfit = funcout(funcindex(ind-1)+1:funcindex(ind));
					yfit = funcout(funcindex(ind-1)+1:funcindex(ind));
				else
%					yfit = funcout(1:funcindex(ind));
					yfit = funcout(1:funcindex(ind));
				end
	
				hold off; errorbar(xdata,ydata,edata,'b--');
				hold on; plot(xfit,yfit,'r-','linewidth',1,[0 80],[0 0],'k--');
%				axis([DAT.xdat(1,ind) DAT.xdat(end,ind)]);
				axis([DAT.eng(1) DAT.eng(end)]);
				title(['column: ',num2str(ind),', q point: ',num2str(qpoint)]);
				xlabel('Energy (meV)');
				ylabel('Intensity (arb. units)');
				legend('Actual Data', 'Fit Data')
				pause
			end
		end
	end
end
hold off

function [F,J,varargout] = objective(SYMS,vars,y_obs);
	disp('Calling objective...')
	weights = 1;		% DerivativeCheck passed when using weights

	if nargout > 1
		[funcout,J]=calc_full_model(SYMS,vars);
	else
		funcout = calc_full_model(SYMS,vars);
	end

	if weights
		F = (funcout - y_obs) .* sqrt(SYMS{1}.VARS.w_obs);
	else
		F = funcout - y_obs;
	end
end
end
