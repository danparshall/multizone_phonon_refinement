function SYMS = validate_user_config(file_CONFIG);
% A user can set the parameters for a refinement within a config file


% energy range to consider
% background file
% which centers to refine (can be overridden in "assemble_VARS" if no Q has decent data)
% which widths to refine (ditto)
% ability to override starting centers and/or widths (e.g. when the DFT calc isn't that good)

% this should apply in "make_AUX_mask" and/or "generate_AUX"


display(['Reading file locations and other config settings...'])
CONFIG = file_CONFIG();

has_preds = length(CONFIG.paths_predictions) > 0;

class_dat = class(CONFIG.paths_data);
switch class_dat

    case 'cell'
    % if paths are arrays, confirm that they are the same length, contents are correct, etc.
        n_dat = length(CONFIG.paths_data);

        for i = 1:n_dat
            assert( isa(CONFIG.paths_data{i}, 'char'), 'All elements of "paths_data" must be strings.');
        end

        if has_preds
            assert(isa(CONFIG.paths_predictions, 'cell'), '"paths_data" and "paths_predictions" must both be strings, arrays of strings');
            assert(length(CONFIG.paths_predictions == n_dat), "You must have a predictions file for every dataset (if using predictions).");
            for i = 1:n_dat
                assert( isa(CONFIG.paths_predictions{i}, 'char'), 'All elements of "paths_predictions" must be strings.'));
            end
        end


    case 'char'
    % if paths are strings, wrap them in cell arrays for convenience
        n_dat = 1;
        tmp = cell();
        tmp{1} = CONFIG.paths_data;
        CONFIG.paths_data = tmp;

        if has_preds
            assert(isa(CONFIG.paths_predictions), 'char'), '"paths_data" and "paths_predictions" must both be strings, arrays of strings');
            tmp = cell();
            tmp{1} = CONFIG.paths_predictions;
            CONFIG.paths_predictions = tmp;
        end

    otherwise
        error('"paths_data" must be a string, or cell array of strings')
end




% Verify that files exist
for i = 1:n_dat
    if ~exist(CONFIG.paths_data{i}, 'file')
        error(["Unable to proceed, DAT file not found at : ", CONFIG.path_data, ""])
    end

    if has_preds && ~exist(CONFIG.paths_predictions{i}, 'file')
        error(["Unable to proceed, DFT predictions file not found at : ", CONFIG.path_predictions])
    end

    if ~exist(CONFIG.paths_background{i}, 'file')
        warning(["Unable to locate background file (but fitting can still proceed) at : ", CONFIG.path_background])
    end
end



% validate lengths of dataset-relevant parameters
assert( (CONFIG.estimate_elastic == 0) || (CONFIG.estimate_elastic == 1), '"estimate_elastic" must be 0 or 1') 
assert(length(CONFIG.e_min) == n_dat), "You must provide a lower bound on fitting energy for each dataset."))
assert(length(CONFIG.e_max) == n_dat), "You must provide an upper bound on fitting energy for each dataset."))


if length(CONFIG.Ei) > 0
    has_Eis = 1;
    assert(length(CONFIG.Ei) == n_dat, "If you specify initial energy, one must be provided for each dataset.")
    assert(length(CONFIG.chopfreq == n_dat), "If you specifiy initial energy, you must provide a chopper frequency for each dataset.")
else
    has_Eis = 0;
end


% Load datasets and predictions
SYMS = load_data_and_preds(CONFIG.path_data, CONFIG.path_predictions, CONFIG.estimate_elastic);
assert(length(SYMS) == n_dat, "ERROR : number of SYMS is not equal to provided number of datasets")


for i_sym = 1:n_dat

    %% verify presence/consistency of initial energy for each dataset
    if isfield(SYMS{i_sym}.DAT, 'Ei');
        Ei_tmp = SYMS{i_sym}.DAT.Ei;
        if has_Eis && (SYMS{i_sym}.DAT.Ei ~= CONFIG.Ei[i_sym])
            display(["Discrepancy between Ei specified in CONFIG file (%f) and DATA file (%f); using value from CONFIG file"], [CONFIG.Ei[i_sym], Ei_tmp]);
            SYMS{i_sym}.DAT.Ei = CONFIG.Ei[i_sym];
        end
    elseif has_Eis
        display(["Datafile doesn't have Ei specified; using value provided in config file"])
        SYMS{i_sym}.DAT.Ei = Ei[i_sym];
    else
        error("Ei not found in DAT file, nor in CONFIG.  Need to specify Ei.")
    end


    %% verify presence/consistency of chopper frequency for each dataset
    if isfield(SYMS{i_sym}.DAT, 'chopfreq');
        if has_chopfreqs && (SYMS{i_sym}.DAT.chopfreq ~= CONFIG.chopfreq[i_sym])
            display(["Discrepancy between chopfreq specified in CONFIG file (%f) and DATA file (%f); using value from CONFIG file"], [CONFIG.chopfreq[i_sym], chopfreq_tmp]);
            SYMS{i_sym}.DAT.chopfreq = CONFIG.chopfreq[i_sym];
        end
    elseif has_Eis
        display(["Datafile doesn't have chopfreq specified; using value provided in config file"])
        SYMS{i_sym}.DAT.chopfreq = chopfreq[i_sym];
    else
        error("chopfreq not found in DAT file, nor in CONFIG.  Need to specify chopfreq.")
    end

end


%% validate number of centers, widths, etc
start_cens = SYMS{1}.startvars(:, 1, 1);
if CONFIG.estimate_elastic
    start_cens = start_cens[2:]
end
n_ph = length(start_cens);

assert(length(CONFIG.refine_centers) == length(CONFIG.refine_widths), '"refine_centers" and "refine_widths" must be the same length.')
has_refine_flags = length(CONFIG.refine_centers) > 0;

if has_refine_flags
    assert(n_ph == length(CONFIG.refine_centers), 'Flags for "refine_centers" must be same length as number of centers in predictions'));
end




