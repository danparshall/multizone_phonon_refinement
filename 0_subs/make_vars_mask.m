function [vars_mask, jac_nnz] = make_vars_mask(AUX, cols_prev);
% For every variable located within auxvars, this function determines the column number it will be assigned in the full jacobian.
% the vars_mask has the same shape as auxvars/freevars, but the nonzero elements are occupied by the column of the full jacobian.
%
% So, something like full_jac(FULL_ROWS, vars_mask(find(freevars))) = aux_jac(AUX_ROWS, find(freevars)) will populate the full jacobian (for this SYM)
% this should also be invertible : auxvars(find(freevars)) = VARS.allvars(vars_mask(find(freevars))) should update auxvars with the values in VARS.allvars

mask = AUX.mask;
freevars = AUX.freevars;
n_cen = size(freevars, 1) - 1;
nQ = size(freevars, 2) - 1;

% for every variable, we store the column index it will have in the jacobian of all variables
vars_mask = zeros(n_cen+1, nQ+1, 2);  % same shape as auxvars

% centers and widths are always the first and second batches of columns, regardless of which SYM this is
cols_cen = find(freevars(1:end-1, 1, 1));             % we check freevars, in case some phonons aren't able to be fit at all
cols_wid = n_cen + find(freevars(1:end-1, 1, 2));     % check explicitly which (if any) peaks are having their widths fit as well
i_cens = [1:n_cen]' .* freevars(1:end-1, 1, 1);
i_cens = i_cens(find(i_cens));                        % retains just indices of phonons we are fitting
i_wids = [1:n_cen]' .* freevars(1:end-1, 1, 2);
i_wids = i_wids(find(i_wids));
vars_mask(i_cens, 1, 1) = cols_cen;
vars_mask(i_wids, 1, 2) = cols_wid;

jac_nnz = length(find(mask)) * (length(i_cens) + length(i_wids)) ;

for i_q = 1:nQ
    cols_hts = cols_prev + (i_q - 1)*(n_cen + 1) + find(freevars(:, 1+i_q, 1));                     % heights and constant BG
    cols_res = cols_prev + (i_q - 1)*(n_cen + 1) + nQ*(n_cen + 1) + find(freevars(:, 1+i_q, 2));    % linear BG (and res widths)

    vars_mask(find(freevars(:, 1+i_q, 1)), 1+i_q, 1) = cols_hts;
    vars_mask(find(freevars(:, 1+i_q, 2)), 1+i_q, 2) = cols_res;

    num_added = length(find(mask(:,i_q))) * (length(cols_hts) + length(cols_res));
    jac_nnz = jac_nnz + num_added;
end  % end Q loop
disp(["est. vars_nnz : ", num2str(jac_nnz)]);
var_check = find(freevars) == find(vars_mask);
assert(sum(var_check) == length(var_check), "vars_mask should only be occupied where freevars==1")
