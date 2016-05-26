function SYM = refine_phonons_multizones(SYM);
% SYM = refine_phonons_multizones(SYM)
%	SYM is an ensemble of data from a single reduced-q point.  
%	
%	Each element of SYM will correspond to data from a different .sqw file,
%	where each .sqw has a different experimental condition (such as different 
%	crystal orientations, incident energies, etc.).
%
%	Each element of SYM contains a data bundle (DAT), as well as auxilary
%	information about that data (AUX).  AUX contains the variables for the
%	model, such as peak centers/heights/wids, instrument resolution, etc.
%	SYM elements also contain Ei and chopper frequency.
%	
%	The first element of SYM also contains VARS, which is an ensemble of AUX
%
%	DAT should have fields for x_dat, y_dat, e_dat (energy, intensity, error), 
%	as well as HKL_vals (and ideally Q_mags).

%% ## !!! may have to install and load 'optim' package:
%	pkg install -forge optim
%	pkg load optim
%	pkg describe -verbose optim		# lists all functions in optim

debug = 1;

VARS=SYM{1}.VARS;
varsin=VARS.varsin(:);
ydatin=VARS.ydatin(:);
wdatin=VARS.wdatin(:);


%display number of Q-points used
nZones = 0;
for ind = 1:length(SYM)
	nZones = nZones + length(find(sum(SYM{ind}.AUX.mask)));
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
									@(vars) objective(SYM,vars,ydatin), ...
									varsin, bounds_lo, bounds_hi, opts);


%% === post-process ===

	report_exitflag(exitflag);

	%	eliminates centers that had no data to fit
%	varsout(find(SYM{1}.VARS.freevars([1:end-1],1,1) == 0)) = NaN;


%% === update and uncertainty ===
SYM=update_AUX(SYM,varsout);
SYM=make_VARS(SYM);

%% check for NaN
%nansum_multizone = sum(isnan(SYM{1}.AUX.auxvars(:,1,1)))


if 0
	unc = calc_unc(SYM);
end



% === plotting and debugging ===
if 0
	for sfile = 1:length(SYM)
		cen = SYM{1}.VARS.allvars(1:end-1,1)
		funcindex = cumsum(sum(SYM{sfile}.AUX.mask));
		mask = SYM{sfile}.AUX.mask;
		SYM = SYM{sfile}.DAT;
		for ind = 1:SYM{sfile}.AUX.Nq
			if sum(mask(:,ind))>0
				qpoint = SYM.HKLvals(ind,:);
				xdata = SYM.xdat(mask(:,ind),ind);
				ydata = SYM.ydat(mask(:,ind),ind);
				edata = SYM.edat(mask(:,ind),ind);
%save_data(xdata',ydata',edata',qpoint);
%SYM{1}.VARS.allvars(:,ind+1,:)
%SYM{1}.AUX.freevars(:,ind+1,:)
				xfit = SYM.xdat(mask(:,ind),ind);
				if ind ~= 1
					yfit = funcout(funcindex(ind-1)+1:funcindex(ind));
				else
					yfit = funcout(1:funcindex(ind));
				end
	
				hold off; errorbar(xdata,ydata,edata,'b--');
				hold on; plot(xfit,yfit,'r-','linewidth',1,[0 80],[0 0],'k--');
				axis([SYM.xdat(1,ind) SYM.xdat(end,ind)]);
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


%% === objective function used for lsqnonlin ===
function [F,J,varargout] = objective(SYM,vars,ydatin);

	weights = 1;		% DerivativeCheck passed when using weights

	if nargout > 1
		[funcout,jacobian]=calc_model_multiQ(SYM,vars);

		if weights
			J = jacobian .* repmat(sqrt(SYM{1}.VARS.wdatin), 1, size(jacobian,2));
		end
	else
		funcout = calc_model_multiQ(SYM,vars);
	end

	if weights
		F = (funcout - ydatin) .* sqrt(SYM{1}.VARS.wdatin);
	else
		F = funcout - ydatin;
	end
end
end
