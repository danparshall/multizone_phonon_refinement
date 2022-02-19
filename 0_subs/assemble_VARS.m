function SYMS=assemble_VARS(SYMS);
% SYMS = assemble_VARS(SYMS);
% Assembles individual AUX.auxvars into VARS.allvars
%
% auxvars are :
%	(:,:,1)=[cen, height;0, constant BG]
%	(:,:,2)=[wid, res; 0, constant BG];

% allvars is a column vector; the first Nph entries are centers;
% the next Nph entries are widths.  After that are the Q-specific
% variables for each SYM.  First (Nph heights + 1 const. BG) for 
% each Q-point, followed by (Nph reswidths + 1 linear BG) for
% each Q-point of the first SYM; subsequent SYMs are appended in
% the same pattern (Nph + 1: for heights/const, then another Nph+1
% for reswidth/linear).  These correspond to the columns of the 
% full Jacobian.
% 
% VARS.indfree contains the indices for those variables which
% are actually fitted; running VARS.allvars(VARS.indfree) will
% return just the subset which is needed for fitting.  Likewise,
% VARS.indfree can be used to select the columns of the jacobian
% which correspond to those variables.
%
% y_obs and w_obs are column vectors which stack the data from
% each SYM.


debug = 0;

n_sym = length(SYMS);
n_cen = size(SYMS{1}.AUX.freevars, 1) - 1;
SYMS{1}.VARS.Nph = n_cen;

%%%% Here we check across all SYMS: in order to fit centers/widths, the peak has to be fittable SOMEWHERE
poss_cen_fitting = zeros(n_cen, 1);
for i_sym = 1:n_sym
    this_poss = sum(SYMS{i_sym}.AUX.freevars(1:n_cen, 2:end, 1), 2) > 0;
    poss_cen_fitting = poss_cen_fitting + this_poss;
end  % end cen-fitting loop
poss_cen_fitting = poss_cen_fitting > 0;

for i_sym = 1:n_sym
    if (length(find(poss_cen_fitting)) ~= n_cen)
        if  (i_sym == 1);
            disp("WARNING : some phonon peaks not fittable in any zone; these are being forced to non-free");
            if system_octave; fflush(stdout); end;
        end;
        SYMS{i_sym}.AUX.freevars(1:n_cen, 1, 1) = poss_cen_fitting .* SYMS{i_sym}.AUX.freevars(1:n_cen, 1, 1);
        SYMS{i_sym}.AUX.freevars(1:n_cen, 1, 2) = poss_cen_fitting .* SYMS{i_sym}.AUX.freevars(1:n_cen, 1, 2);
        if debug; disp(SYMS{i_sym}.AUX.freevars); end;
    end 

    %% make mapping for jacobian indices.  (we do it here rather than "generate_AUX", since we may be changing freevars above)
    % mapping from fitted variables to the indices corresponding to those variables within AUX.auxjac (fulljac calculated elsewhere)
    shape_fvars = size(SYMS{i_sym}.AUX.freevars);
    SYMS{i_sym}.AUX.jac_inds = reshape([1:prod(shape_fvars)], shape_fvars) .* SYMS{i_sym}.AUX.freevars;

    % mapping from energy points contributing to the fitting, to their indices
    shape_eng = size(SYMS{i_sym}.AUX.mask);
    SYMS{i_sym}.AUX.eng_inds = reshape([1:prod(shape_eng)], shape_eng) .* SYMS{i_sym}.AUX.mask;
end


%%%% Given our updated SYMS, create the composite VARS

nQs = [];   % number of Q-points for which each SYM has data
nEs = [];   % length of energy array for each SYM
y_obs = [];
w_obs = [];
func_mask = [];
jac_nnz = [];
alt_nnz = [];

obs_prev = 0;
cols_prev = 2*n_cen;			% the first 2*n_cen cols of the full jacobian are reserved for cens & wids
for i_sym = 1:n_sym;
    if debug; disp(''); disp(["Assembling VARS using data from SYM : " num2str(i_sym)]); end;

    AUX = SYMS{i_sym}.AUX;
    DAT = SYMS{i_sym}.DAT;
    nQ = AUX.Nq;
    nE = length(AUX.eng);

    y_masked = DAT.y_dat .* AUX.mask;
    y_obs = [y_obs; y_masked(:)];

    w_masked = (1./DAT.e_dat).^2 .* AUX.mask;
    w_obs = [w_obs; w_masked(:)];

    this_mask = find(AUX.mask);
    func_mask = [func_mask; obs_prev + this_mask(:)];
    obs_prev = obs_prev + numel(AUX.mask);

    % append number of energy/Q points to our running total
    % nQ x nE is the nominal number of rows from this SYM (i.e. nominal number of rows in the Jacobian)
    nQs = [nQs nQ];
    nEs = [nEs nE];

    if i_sym == 1
        if debug; disp('initializing VARS...'); end;
        % the first 2*Nph elements of VARS are cens & wids
        bounds_L = [AUX.bounds_L(1:end-1, 1, :)]; bounds_L = bounds_L(:);
        bounds_H = [AUX.bounds_H(1:end-1, 1, :)]; bounds_H = bounds_H(:);
        freevars = [AUX.freevars(1:end-1, 1, :)]; freevars = freevars(:);
        allvars = [AUX.auxvars(1:end-1, 1, :)];  allvars = allvars(:);
        prev_vars = length(freevars);
    end  % end if


    % for every element of freevars, identify the corresponding column index of the full jacobian
    [vars_mask, aux_nnz] = make_vars_mask(AUX, cols_prev);
    SYMS{i_sym}.AUX.vars_mask = vars_mask;
    SYMS{i_sym}.AUX.aux_nnz = aux_nnz;
	cols_added = 2 * AUX.Nq * (AUX.Nph + 1);
    cols_prev = cols_prev + cols_added;


    % each SYM appends Nq*(Nph+1) height+constant, and the same number of reswid+linear
    tmp = AUX.bounds_L(:, 2:end, :);
    bounds_L = [bounds_L; tmp(:)];
    tmp = AUX.bounds_H(:, 2:end, :);
    bounds_H = [bounds_H; tmp(:)];
    tmp= AUX.freevars(:, 2:end, :);
    freevars = [freevars; tmp(:)];
    tmp= AUX.auxvars(:, 2:end, :);
    allvars = [allvars; tmp(:)];
    prev_vars = [prev_vars, length(freevars)];
    jac_nnz = [jac_nnz, SYMS{i_sym}.AUX.aux_nnz];


    % alternate method to determine NNZ per SYM
    if debug
        alt = 0;
        for iq = 1:nQ
            free_eng = AUX.mask(:,iq);
            free_vars = [reshape(AUX.freevars(:, 1, :), [], 1); reshape(AUX.freevars(:, 1+iq, :), [], 1)]';  % column vector
            alt = alt + length(find(free_eng * free_vars));
        end
        alt_nnz = [alt_nnz, alt];
    end % end alt-NNZ calculator
end  % end first SYM loop

assert(sum(nQs .* nEs) == length(y_obs), "The output function seems to be the wrong length.")

if debug
    disp(["VARS.jac_nnz : ", num2str(jac_nnz)])
    disp(["VARS.alt_nnz : ", num2str(alt_nnz)])
    disp(["VARS freevars has length " num2str(length(freevars)) ", with " num2str(length(find(freevars))) " fittable"]);
end

SYMS{1}.VARS.allvars = allvars;
SYMS{1}.VARS.bounds_L = bounds_L;
SYMS{1}.VARS.bounds_H = bounds_H;
SYMS{1}.VARS.freevars = freevars;
SYMS{1}.VARS.indfree = find(freevars);
SYMS{1}.VARS.varsin = SYMS{1}.VARS.allvars(SYMS{1}.VARS.indfree);
SYMS{1}.VARS.nQs = nQs;
SYMS{1}.VARS.nEs = nEs;
SYMS{1}.VARS.y_obs = y_obs;
SYMS{1}.VARS.w_obs = w_obs;
SYMS{1}.VARS.func_mask = func_mask;
SYMS{1}.VARS.prev_vars = prev_vars;
SYMS{1}.VARS.jac_nnz = jac_nnz;
if debug; SYMS{1}.VARS.alt_nnz = alt_nnz; end;
