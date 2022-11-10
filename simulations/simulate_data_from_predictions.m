function [SYM,startvars,true_y] = simulate_data_from_predictions(SYM, sim_vars, junk_scale, pred_error)



function [modelout, jacout] = sim_escan(eng, cens, heights, wid_ph, wid_res, asymm)
    eng=eng(:)';
    modelout=zeros(size(eng));

    for ind=1:N_ph;
        [splitgauss,splitjac] = calc_splitgauss_full(eng, cens(ind), heights(ind), wid_ph(ind), wid_res(ind), asymm);
        modelout = modelout + splitgauss';
    end
end


%% 
eng = SYM.eng;
N_q = size(sim_vars, 2) - 2
N_ph = size(sim_vars, 1)
cens = sim_vars(:,1);

%% resolution-adjusted widths
reswids=merchop(SYM.Ei, SYM.chopfreq, cens);
wids_ph = sim_vars(:,2);
asymm = 1.7;


heights = sim_vars(:, 3:end);


if ~exist('junk_scale'); junk_scale=0.4; end
sim_y = zeros(length(eng),N_q);
sim_e = sim_y;
true_y = sim_y;
sim_x = repmat(eng,1,size(sim_y,2));	% simulating Q

rand('seed', 42);
for inq=1:N_q

	% model of ideal data
    hts = heights(:, inq);   % heights at just this Q
	[model] = sim_escan(eng, cens, hts, wids_ph, reswids, asymm);
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
	xdat(xdat > SYM.E_maxes(inq)) = 0;
	ydat(xdat == 0) = 0;
	edat(xdat == 0) = 0;

	% update
	sim_x(:,inq)=xdat;
	sim_y(:,inq)=ydat;
	sim_e(:,inq)=edat;
end

SYM.DAT.y_dat = sim_y;
SYM.DAT.e_dat = sim_e;
SYM.DAT.x_dat = sim_x;
SYM.DAT.eng = eng;


if 1
	% make starting variables by adding noise to true variables (since our DFT prediction is never perfect)
	if ~exist('pred_error')
		disp('No pred_error specified; using junk_scale')
		pred_error = junk_scale/4;
	end

	% 95% of the time, prediction starts within $PRED_ERROR percent of true value
	cen_noise = cens .*  (pred_error*randn(N_ph, 1));
	wid_noise = sim_vars(:,2) .* (pred_error*randn(N_ph, 1));
	hts_noise = heights .* (pred_error*randn(N_ph, N_q));

else
	% don't add noise; simulated DFT predictions correspond exactly to the true underlying values
	cen_noise = zeros(N_ph, 1);
	wid_noise = zeros(N_ph, 1);
	hts_noise = zeros(N_ph, N_q);
end  % end noise addition

start_cens = sim_vars(:,1) + cen_noise;
start_wids = sim_vars(:,2) + wid_noise;
start_hts = heights + hts_noise;
startvars = [start_cens, start_wids, start_hts ];
end