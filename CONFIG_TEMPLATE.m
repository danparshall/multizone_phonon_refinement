function CONFIG = CONFIG;
%	Template file for user preferences for refining individual symmetry bundles.
%	Similar to an EXP file used for ResLib or SNAXS
%
%   Since MPR can refine multiple symmetry bundles simultaneously, the parameters here
%   need to agree.  E.g., the number of elements in "paths_data" must be equal to that
%   in "paths_predictions", the length of "refine_centers" must be equal to that of 
%   "refine_widths"


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% DATASETS
%%% parameters related to datasets. Number of elements of all of these must be equal.

% filepath for DAT file.  Can be string, or cellarray of strings.
CONFIG.paths_data = "~/path/to/data/file";


% filepath for predictions file.  Can be string, or cellarray of strings.
CONFIG.paths_predictions = "~/path/to/predictions/file";


% filepath for background file.  Can be string, or cellarray of strings.
CONFIG.paths_background = "~/path/to/background/file";


% whether or not to include estimate of elastic line
CONFIG.estimate_elastic = 1;

% minimum energy to fit
CONFIG.e_min = 5;

% maximum energy
CONFIG.e_max = 70;


% Ei and chopper freq are optional; they will override the values recorded in DAT file if present
CONFIG.Ei = 100;
CONFIG.chopfreq = 600;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% MATERIAL
%%% parameters related to material. Number of elements of all of these must be equal.

% boolean indicating if the center of each phonon should be refined or not
% if empty, MPR will attempt to refine all centers
CONFIG.refine_centers = [0 1 1 1 1 1 0];

% boolean indicating if the width of each phonon should be refined or not
% if empty, MPR will attempt to refine all widths
CONFIG.refine_widths = [0 1 1 1 1 1 0];

% these can override the values from the predictions file; if not empty must be same length as refine_centers
CONFIG.starting_centers = nan;

% ditto
CONFIG.starting_widths = nan;

