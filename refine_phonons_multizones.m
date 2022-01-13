function SYMS = refine_phonons_multizones(SYMS);
% SYMS = refine_phonons_multizones(SYMS)
%	SYMS is an ensemble of data from a single reduced-q point.  
%	
%	Each element of SYMS will correspond to data from a different .sqw file,
%	where each .sqw has a different experimental condition (such as different 
%	crystal orientations, incident energies, etc.).
%
%	Each element of SYMS contains a data bundle (DAT), as well as auxilary
%	information about that data (AUX).  AUX contains the variables for the
%	model, such as peak centers/heights/wids, instrument resolution, etc.
%	SYMS elements also contain Ei and chopper frequency.
%	
%	The first element of SYMS also contains VARS, which is an ensemble of AUX
%
%	DAT should have fields for x_dat, y_dat, e_dat (energy, intensity, error), 
%	as well as HKL_vals (and ideally Q_mags).

%% ## !!! may have to install and load 'optim' package:
%	pkg install -forge optim
%	pkg load optim
%	pkg describe -verbose optim		# lists all functions in optim

debug = 1;

VARS=SYMS{1}.VARS;
varsin=VARS.varsin(:);
ydatin=VARS.ydatin(:);
wdatin=VARS.wdatin(:);


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
		'Jacobian', 'on', ...
		'Display', 'iter');

	
	%% Check that array sizes 
	[~,maxsize] = computer;
	if system_octave;
		% https://github.com/calaba/octave-3.8.2-enable-64-ubuntu-14.04
		if sizemax < 2^31;
			warning(' Octave should probably be recompiled to handle larger arrays, see https://github.com/calaba/octave-3.8.2-enable-64-ubuntu-14.04');
		end
	end
	assert(length(varsin)*length(ydatin) < maxsize, 'The number of elements in the Jacobian is greater than the maxiumum array length that can be indexed.');


	if debug
		disp([' size(vars):' num2str(size(varsin))]);
		disp([' size(ydat):' num2str(size(ydatin))]);
		disp('  OPTS : ')
		opts
	end


	bounds_lo = VARS.bndsLO(VARS.indfree);
	bounds_hi = VARS.bndsHI(VARS.indfree);
[varsout,resnorm,resid,exitflag] = lsqnonlin( ...
									@(vars) objective(SYMS,vars,ydatin), ...
									varsin, bounds_lo, bounds_hi, opts);


%% === post-process ===

	report_exitflag(exitflag);

	%	eliminates centers that had no data to fit
%	varsout(find(SYMS{1}.VARS.freevars([1:end-1],1,1) == 0)) = NaN;


%% === update and uncertainty ===
SYMS=update_AUX(SYMS,varsout);
SYMS=make_VARS(SYMS);

%% check for NaN
%nansum_multizone = sum(isnan(SYMS{1}.AUX.auxvars(:,1,1)))


if 0
	unc = calc_unc(SYMS);
end


%% === below here is just for displaying graphs and debugging

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
		cen = SYMS{1}.VARS.allvars(1:end-1,1)
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

function [F,J,varargout] = objective(SYMS,vars,ydatin);

	if nargout > 1
		[funcout,J]=calc_model_multiQ(SYMS,vars);
	else
		funcout = calc_model_multiQ(SYMS,vars);
	end
	F = funcout - ydatin;
end
end
