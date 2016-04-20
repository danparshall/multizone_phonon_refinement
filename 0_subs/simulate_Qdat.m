function [simy,sime,htout,model,junk]=simulate_Qdat(gvars,Nsim,eng,junkscale);
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
%	model is the "true" intensity data, without noise
%	junk is the actual noise

if ~exist('junkscale'); junkscale=0.4; end
simy=zeros(length(eng),Nsim);
sime=simy;

htout=zeros(size(gvars,1),Nsim);

for inq=1:Nsim

	% gvars is center, height, Lwidth, Rwidth for each peak
	% this sets heights
	rsim=rand([1 size(gvars,1)]);
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

