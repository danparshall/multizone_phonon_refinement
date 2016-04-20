function SYMDAT = refine_phonons_multizones(SYMDAT);
% SYMDAT = refine_phonons_multizones(SYMDAT)
%	SYMDAT is an ensemble of data from a single reduced-q point.  
%	
%	Each element of SYMDAT will correspond to data from a different .sqw file,
%	where each .sqw has a different experimental condition (such as different 
%	crystal orientations, incident energies, etc.).
%
%	Each element of SYMDAT contains a data bundle (SQW), as well as auxilary
%	information about that data (AUX).  AUX contains the variables for the
%	model, such as peak centers/heights/wids, instrument resolution, etc.
%	SYMDAT elements also contain Ei and chopper frequency.
%	
%	The first element of SYMDAT also contains VARS, which is an ensemble of AUX
%
%	SQW should have fields for x_dat, y_dat, e_dat (energy, intensity, error), 
%	as well as HKL_vals (and ideally Q_mags).

% Perhaps should retitle these things:
%	SQW 	-> DAT
%	VARS 	-> AUXS
%	SYMDAT	-> SYM



VARS=SYMDAT{1}.VARS;
varsin=VARS.varsin;
ydatin=VARS.ydatin;
wdatin=VARS.wdatin;

options.bounds=[VARS.bndsLO(VARS.indfree) VARS.bndsHI(VARS.indfree)];

stol=0.00001;
niter=10000;
dp=0.00001*ones(size(varsin));
%dFdp='calc_JAC_multiQ';
varsin=varsin(:);
ydatin=ydatin(:);
wdatin=wdatin(:);

const = SYMDAT{1}.VARS.allvars(end,2:end,1);
lin = SYMDAT{1}.VARS.allvars(end,2:end,2);

%display number of Brillouin zones used
nZones = 0;
for ind = 1:length(SYMDAT)
	nZones = nZones + length(find(sum(SYMDAT{ind}.AUX.mask)));
end
disp([' Used ' num2str(nZones) ' seperate Brillouin zones']);


%% === fit data ===
[funcout,varsout,cvg,iter,corp,covp,covr,stdres,Z,r2]=leasqr(SYMDAT,ydatin,varsin,...
		'calc_DAT_multiQ',stol,niter,wdatin,dp,'calc_JAC_multiQ',options);


%[funcout,varsout,cvg,iter,corp,covp,covr,stdres]=leasqr(SYMDAT,ydatin,varsin,...
%		'calc_DAT_multiQ',stol,niter,wdatin,dp);



%% === post-process ===
	if cvg~=1;
		disp(' WARNING in "refine_phonons_multizones" : no covergence');
	end

	%	eliminates centers that had no data to fit
	varsout(find(SYMDAT{1}.VARS.freevars([1:end-1],1,1) == 0)) = NaN;


%% === update and uncertainty ===
SYMDAT=update_AUX(SYMDAT,varsout);
unc = calc_unc(SYMDAT);



%% === below here is just for displaying graphs and debugging

% === plotting ===
if 0
	for sfile = 1:length(SYMDAT)
		cen = SYMDAT{1}.VARS.allvars(1:end-1,1)
		funcindex = cumsum(sum(SYMDAT{sfile}.AUX.mask));
		mask = SYMDAT{sfile}.AUX.mask;
		SYM = SYMDAT{sfile}.SQW;
		for ind = 1:SYMDAT{sfile}.AUX.Nq
			if sum(mask(:,ind))>0
				qpoint = SYM.HKLvals(ind,:);
				xdata = SYM.xdat(mask(:,ind),ind);
				ydata = SYM.ydat(mask(:,ind),ind);
				edata = SYM.edat(mask(:,ind),ind);
%save_data(xdata',ydata',edata',qpoint);
%SYMDAT{1}.VARS.allvars(:,ind+1,:)
%SYMDAT{1}.AUX.freevars(:,ind+1,:)
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
