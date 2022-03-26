function SYMS = load_data_and_preds(paths_DAT, paths_PRED, estimate_elastic, varargin);

if ~exist('estimate_elastic')
    estimate_elastic = 1;
end

class_dat = class(paths_DAT);
class_pred = class(paths_PRED);

switch class_dat
    case 'char'
        n_dat = 1;
        tmp = cell();
        tmp{1} = paths_DAT;
        paths_DAT = tmp;
    case 'cell'
        n_dat = length(paths_DAT);
        for i = 1:n_dat
            assert( isa(paths_DAT{i}, 'char'), 'All elements of "paths_DAT" must be strings.')
        end
    otherwise
        error('ERROR : "paths_DAT" must be a string, or cell array of strings')
end


switch class_pred
    case 'char'
        n_pred = 1;
        tmp = cell();
        tmp{1} = paths_PRED;
        paths_PRED = tmp;
    case 'cell'
        n_pred = length(paths_PRED);
        for i = 1:n_pred
            assert( isa(paths_PRED{i}, 'char'), 'All elements of "paths_PRED" must be strings.')
        end
    otherwise
        error('ERROR : "paths_PRED" must be a string, or cell array of strings')
end

SYMS = cell();
assert(n_dat == n_pred, "Need to have a prediction file for every DAT file.")
for i_file = 1:n_dat
    DAT = load_DAT_file(paths_DAT{i_file});
    startvars = load_pred_file(paths_PRED{i_file}, DAT);

%    if exist('censcale')
%        disp(["Scaling phonon energies by ", num2str(censcale)])
%        startvars(:, 1) = censcale * startvars(:, 1);
%    end

    if estimate_elastic
        disp('Adding estimate of elastic line to fitting');
        elastic = estimate_elastic_line(DAT);
        startvars = [ [0, 1, elastic(:)']; startvars];
    end
    SYMS{i_file}.DAT = DAT;
    SYMS{i_file}.startvars = startvars;
end
