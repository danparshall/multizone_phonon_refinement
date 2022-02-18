function startvars = load_pred_file(path_to_pred_file, DAT, varargin);
% Reads a file containing predictions for phonon centers, widths, and heights.
% Conducts basic consistency checking, so even if a file is hand-edited it
% won't cause problems downstream.
%
% TODO: support for spreadsheets (which could store Q value as column name)

is_textfile = 0;

if exist(path_to_pred_file, 'file') == 2
    [fPath, fName, fExt] = fileparts(path_to_pred_file);
    switch lower(fExt)
        case '.mat'
            % matlab-compatible file
            PRED = load(path_to_pred_file);

        case '.txt'
            % text file.  Can't store metadata, so assumed to be a CSV in proper format: [cens, wids, heights]
            is_textfile = 1;

            tmp = csvread(path_to_pred_file);
            PRED.centers = tmp(:,1);
            PRED.widths = tmp(:,2);
            PRED.heights = tmp(:, 3:end);

            % make dummy Q_hkl
            n_q = size(tmp, 2) - 2;
            qs = [1:n_q]';
            PRED.Q_hkl = [zeros(length(qs), 2), qs];

        case '.ods'
            % Open Document Spreadsheet (can have column headings, but only openable with Matlab)
            error("Support for spreadsheets not yet implemented.  Ensure columns are correct, save as .csv, try again.")
            
        otherwise
            error('Unexpected file extension: %s', fExt);
    end
else
    error('File does not exist : %s', path_to_pred_file)
end

disp(" Loading predictions ...")
disp(["    ... " num2str(length(PRED.Q_hkl)) " Q-points"])
disp(["    ... " num2str(length(PRED.centers)) " phonon centers"])


%% == check that file is internally consistent, i.e. that we have predictions for every phonon at every Q-point
%assert( mod(length(PRED.centers), 3) == 0, "There should be 3*N_atom phonons at every Q")      % there ARE, but sometimes degenerate
assert(length(PRED.centers) == length(PRED.widths), "Each phonon must have a corresponding width.")
assert(length(PRED.centers) == size(PRED.heights, 1), "Need height prediction for every phonon")
assert(length(PRED.Q_hkl) == size(PRED.heights, 2), "Must have height values for every Q")


 %% == if we've passed in a DAT struct, check for consistency there as well
if nargin > 1
    disp(['Validating params found in prediction file to those in the DAT struct...'])
    if ~isa(DAT, 'struct')
        disp(["ERROR : DAT is not a struct object"])
    end

    % confirm that n_Q are the same
    assert(size(DAT.y_dat, 2) == length(PRED.Q_hkl), ['    ... number of Q-points not the same (DAT: ' num2str(size(DAT.y_dat, 2)) ', PRED: ' num2str(length(PRED.Q_hkl)) ")"])

    % confirm that actual Q_hkl values are the same
    if is_textfile
        % since the text file can't store metadata, just use the values from the DAT file
        PRED.Q_hkl = DAT.Q_hkl;
    end
    assert(DAT.Q_hkl == PRED.Q_hkl, '    ... the values of Q are not the same');

    % report that everything is okay
    disp(['    ... params are consistent.'])


    ENG_MINIMUM = 5;
    disp([' Adjusting scale of predictions to match data.  Using energies > ' num2str(ENG_MINIMUM)])

    ind_eng = find(DAT.eng > ENG_MINIMUM);
    max_dat = max(DAT.y_dat(ind_eng, :));   % max val in each column
    max_pred = max(PRED.heights);
    ratio = mean(max_dat ./ max_pred);
    PRED.heights = PRED.heights * ratio;
end

%% == assemble 
startvars = [PRED.centers(:), PRED.widths(:), PRED.heights];