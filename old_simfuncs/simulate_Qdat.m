function [simy,sime,htout,startvars]=simulate_Qdat(gvars,Nsim,eng,junkscale);
% Simulate noisy data.  Models idealized data, then adds junk/noise.
%
% INPUTS:
% 	gvars is center, height, Lwidth, Rwidth for each peak
% 	Nsim is number of Q to simulate
% 	eng is energy (x-axis data)
% 	junkscale is scale factor for added noise
%
% OUTPUTS:
%	simy is simulated intensity data
%	sime is errorbars for the same
%	htout is the "true" height of the peaks, without noise
%	startvars : initial fitting params, such as might be generated from SNAXS

%	model is the "true" intensity data, without noise
%	junk is the actual noise

if ~exist('junkscale'); junkscale=0.4; end
simy=zeros(length(eng),Nsim);
sime=simy;

htout=zeros(size(gvars,1),Nsim);
Nph = size(gvars, 1)

for inq=1:Nsim

	% gvars is center, height, Lwidth, Rwidth for each peak
	% this generates random heights for each phonon at this Q-point
	rsim=rand([Nph, 1]);
	gvars(:,2)=rsim/max(rsim);

	% model of ideal data
	[model]=model_multigauss(gvars,eng);

	% generate noise
	junk=randn(size(model));
	junk=junk*junkscale;

	% add noise
	ydat=model+junk;
	edat=sqrt(junk.^2);
	edat=mean(edat)*ones(size(edat));

	
	simy(:,inq)=ydat;
	sime(:,inq)=edat;
	htout(:,inq)=gvars(:,2);
end

% now we make starting values - normally this would come from SNAXS or other calculation
cen_noise = gvars(:, 1) + 0.05 * randn([Nph, 1]);
ht_noise = htout + (0.1 * randn(size(htout)));
wd_noise = 0.01 * abs(randn([Nph, 1])) ;
startvars = [cen_noise, wd_noise, ht_noise];

