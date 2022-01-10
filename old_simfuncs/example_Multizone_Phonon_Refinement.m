function [fitcens,tocs,SNR]=example_Multizone_Phonon_Refinement(split,resfactor,junkscale,N_q22)
% This function provides a demonstration of Multizone Phonon Refinement.
% Simulated data are fed into the MPR solver.  Various system/experimental 
% parameters can be set in order to evaluate the accuracy of the fitting.
%
% MPR assumes that all the data are from the same symmetry point, and thus all
% the phonons have the same energy and intrinsic linewidth.
% 
% It is possible to collect multiple data sets of the same crystal (either in
% different orientations, or using different Ei). Each data set is considered a
% single instance of S(q,w).  MPR can co-refine multiple DAT.
%
% Simulated data is a pair of overlapping peaks with varying intensity in 
% multiple Brillouin zones.  
%
% INPUTS:
%	split 		- separation between the two overlapping peaks
%	resfactor 	- changes resolution width
%	junkscale 	- sets signal/noise ratio
%	N_q22 		- number of Q for which data have been collected

tic
if nargin < 4
	N_q22=5;
	if nargin < 3;
		junkscale=0.1;
		if nargin < 2
			resfactor=1;
			if nargin <1
				split = .1;
			end
		end
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATED DATA

N_q22=9;
N_q52=9;

gvars=[];

center=13;
cens=[center-(split/2) center+(split/2)]';

%cens = [12.95 13.05 22.95 23.05]'


if 1
	Ei=21.9;
	chopfreq=300;
%	eng1=[0:.25:sum(cens)];eng1=eng1(:);
	eng1=[11:.25:15];eng1=eng1(:);
	reswids=merchop(Ei,chopfreq,cens);	reswids=reswids(:)*resfactor;

	% gvars is center, height, Lwidth, Rwidth for each peak
	gvars=[cens zeros(size(cens)) 0.66*reswids 0.33*reswids];

	[ydat,edat,htout,model,junk]=simulate_Qdat(gvars,N_q22,eng1,junkscale);
	xdat=repmat(eng1,1,size(ydat,2));	% simulating Q

	DAT.xdat=xdat;
	DAT.ydat=ydat;
	DAT.edat=edat;
	DAT.htout=htout;
	SYM{1}.DAT=DAT;
	SYM{1}.Ei=Ei;
	SYM{1}.chopfreq=chopfreq;
%	save('TESTDAT2','DAT');

	if 1		% turn this on to add a second SYM
		Ei=52;
		chopfreq=300;
		eng2=[10:.5:20];eng2=eng2(:);
		reswids=merchop(Ei,chopfreq,cens);	reswids=reswids(:);
		gvars=[cens zeros(size(cens)) 0.66*reswids 0.33*reswids];

		[ydat,edat,htout,model2,junk2]=simulate_Qdat(gvars,N_q52,eng2,junkscale);
		xdat=repmat(eng2,1,size(ydat,2));	% simulating Q
		DAT.xdat=xdat;
		DAT.ydat=ydat;
		DAT.edat=edat;
		DAT.htout=htout;
		SYM{2}.DAT=DAT;
		SYM{2}.Ei=Ei;
		SYM{2}.chopfreq=chopfreq;
		junk=[junk(:); junk2(:)];
		model=[model(:); model2(:)];
	else
		disp(' Setting SYM{2}=SYM{1}')
		SYM{2}=SYM{1};
	end
else
	DAT=load('TESTDAT2');
	DAT=DAT.DAT;
	SYM{1}.DAT=DAT;
disp(' LOADED...')
end

if system_octave;
	fflush(stdout);
end

startvars=gvars;
SNR=sum(model(:))/sum(abs(junk(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%size(SYM{1}.DAT.xdat)

%% FITTING METHOD

%% generate AUX and VARS structures prior to fitting
%% here startvars is given, but in general it would come from a SNAXS calculation.
SYM=generate_AUX(SYM,startvars);
SYM=make_VARS(SYM);

% fit data
SYM = refine_phonons_multizones(SYM);

fitcens=SYM{1}.VARS.allvars(:,1);


fitcens=fitcens(1:size(gvars,1))
gvars
if system_octave;
	fflush(stdout);
end


errfit=(fitcens-gvars(:,1))./gvars(:,1);
%errfit=errfit(:)'
fitcens=fitcens(:)';

tocs=toc;
SYMmult=SYM;
