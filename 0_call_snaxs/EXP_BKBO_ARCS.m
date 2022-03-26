function EXP=demoEXP_tof
%===============================================================================
% function EXP=demoEXP_tof
% adaption from ResLib. Q-resolution is not considered, so only need to indicate
% the parameters listed in this file.
%
%	Required:
%		experiment_type='tof'
%		efixed
%		infin
%		sample.a/b/c
%
%	Optional:
%		basis_user			% can switch between various bases
%		calculation_path	% if provided, SNAXS will load this automatically
%		instrument
%===============================================================================

% === for neutron time-of-flight, very little information is needed ===
EXP.experiment_type = 'tof';	% if not specified, SNAXS defaults to 'tas'
EXP.efixed=100;					% Fixed neutron energy in meV
EXP.chopfreq = 600;				% Fermi chopper frequency (sets resolution)
EXP.infin=1;					% Fixed energy (pos. for fixed Ef, neg. for fixed Ei)
EXP.instrument='ARCS';			% Used to determine Energy resolution (and someday, Q-resolution). 
								% Only ARCS is currently implemented. Hope to add MERLIN soon.


if 0
	EXP.calculation_path='DFTs/2x2x1_tilt/POSCAR';
	EXP.dim = [1 1 1];
	EXP.basis_user = [0.5 0.5 0; 0.5 -0.5 0; 0 0 0.5];
elseif 0
	EXP.calculation_path='DFTs/2x2x1_tiltrot/POSCAR';
	EXP.dim = [1 1 1];
	EXP.basis_user = [0.5 0.5 0; 0.5 -0.5 0; 0 0 0.5];
else
%	EXP.calculation_path='/home/dep/0_SNAXS/phonopy_files/example/CaTiO3/POSCAR';
    EXP.calculation_path = '/home/dan/0_sync/Research/0_software/phonopy/example/CaTiO3/POSCAR';
end

EXP.sample.a=4.287;		%A: angstroms
EXP.sample.alpha=90;
%EXP.user_basis = [0.5 0.5 0; -0.5 0.5 0; 0 0 1];
