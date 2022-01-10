function [SYM,startvars,true_y] = simulate_data_from_predictions(SYM, sim_vars, junk_scale)


wid_ratio = 0.33;
function [modelout, jacout] = sim_escan(eng, cens, wids, heights, wid_ratio)
    eng=eng(:)';
    modelout=zeros(size(eng));

    hwhm_1 = wids*wid_ratio;
    hwhm_2 = wids*(1-wid_ratio);

    for ind=1:N_ph;
        [splitgauss,splitjac]=calc_splitgauss_JAC_fast(eng, cens(ind), heights(ind), hwhm_1(ind), hwhm_2(ind));
        modelout = modelout+splitgauss';
    end
end


%% 
eng = SYM.eng;
N_q = size(sim_vars, 2) - 2
N_ph = size(sim_vars, 1)
cens = sim_vars(:,1);

%% resolution-adjusted widths
reswids=merchop(SYM.Ei, SYM.chopfreq, cens);
wids = sqrt(sim_vars(:,2).^2 + reswids.^2);

heights = sim_vars(:, 3:end);


if ~exist('junk_scale'); junk_scale=0.4; end
sim_y=zeros(length(eng),N_q);
sim_e=sim_y;
true_y = sim_y;
sim_x = repmat(eng,1,size(sim_y,2));	% simulating Q

rand('seed', 42);
for inq=1:N_q

	% model of ideal data
    hts = heights(:, inq);   % heights at just this Q
	[model] = sim_escan(eng, cens, wids, hts, wid_ratio);
    true_y(:,inq) = model;

	% generate noise
	junk=randn(size(model));
	junk=junk*junk_scale;

	% add noise
	ydat=model+junk;
	edat=sqrt(junk.^2);
	edat=mean(edat)*ones(size(edat));

	% mask out inaccessible energies
	xdat = eng(:)';
%	xdat(xdat > SYM.E_maxes(inq)) = nan;
%	ydat(isnan(xdat)) = nan;
%	edat(isnan(xdat)) = nan;

	% update
	sim_x(:,inq)=xdat;
	sim_y(:,inq)=ydat;
	sim_e(:,inq)=edat;
end

SYM.DAT.ydat = sim_y;
SYM.DAT.edat = sim_e;
SYM.DAT.xdat = sim_x;


% make starting variables by adding noise to true variables (since our DFT prediction is never perfect)
cen_noise = cens .* randn(N_ph, 1) * 0.15;    % 95% of the time, start within 15% of prediction
wid_noise = sim_vars(:,2) .* (1 + randn(N_ph, 1)*0.1);
hts_noise = heights .* (1 + randn(N_ph, N_q)*0.2);
start_cens = sim_vars(:,1) + cen_noise;
start_wids = sim_vars(:,2) + wid_noise;
start_hts = heights + hts_noise;
startvars = [start_cens, start_wids, start_hts ];

end