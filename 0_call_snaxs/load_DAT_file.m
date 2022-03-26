function DAT = load_DAT_file(path_to_DAT_file)
% Load datafile with experimental results; do quick sanity-checking.

if exist(path_to_DAT_file, 'file') == 2
    [fPath, fName, fExt] = fileparts(path_to_DAT_file);
    switch lower(fExt)
        case '.mat'
            % matlab-compatible file
            DAT = load(path_to_DAT_file);
  
        otherwise
            error('Unexpected file extension: %s', fExt);
    end
else
    error('File does not exist : %s', path_to_DAT_file)
end


assert(size(DAT.e_dat) == size(DAT.y_dat))
n_E = size(DAT.y_dat, 1);
n_Q = size(DAT.y_dat, 2);


%% == check the Q-values ===
if isfield(DAT, 'HKL_vals')
    if ~isfield(DAT, 'Q_hkl')
        % some batch of code was using the wrong label briefly.  Correct the field name
        DAT.Q_hkl = DAT.HKL_vals;
        rmfield(DAT, 'HKL_vals');
    else
        % Not sure how there could be both.  But it's probably okay as long as they're identical.
        assert(DAT.Q_hkl == DAT.HKL_vals);
    end
end
assert(n_Q == length(DAT.Q_hkl), "Must have HKL defined for every Q-point.");


%% == check the energy values ==
if isfield(DAT, 'x_dat')

    if size(DAT.x_dat, 2) == 1
        % just a column array
        xlen = length(DAT.x_dat);
    elseif size(DAT.x_dat, 1) == 1
        % okay, we'll MAKE it a column array
        DAT.x_dat = DAT.x_dat(:);
        xlen = length(DAT.x_dat);
    elseif size(DAT.x_dat, 2) == size(DAT.y_dat, 2)
        % energy array for every Q - check if the energies are the same everywhere, raise an error if not
        xlen = size(DAT.x_dat, 1);
        mpr_eng_error = "MPR assumes the energy grid within a single SYM is the same in every Q.  You can update it with a couple days work, or make a separate SYM for each energy grid."
        assert(mean(DAT.x_dat, 2) == DAT.x_dat(:, 1), mpr_eng_error);
    else
        % something has gone wrong
        error("ERROR : size of energy array seems inconsistent");
    end


    if isfield(DAT, 'eng')
        DAT.eng = DAT.eng(:);
        assert(length(DAT.eng) == xlen, "Inconsistent length: xdat and eng");
    else
        DAT.eng = DAT.x_dat(:, 1);
    end

    rmfield(DAT, 'x_dat');
end
assert(isfield(DAT, 'eng'), "DAT struct must have 'eng' field.")
assert(sum(~isfinite(DAT.eng)) == 0, "Problem with energy grid.  All values must be finite.")
assert(length(DAT.eng) == size(DAT.y_dat, 1), "Energy grid is inconsistent with y data.")

