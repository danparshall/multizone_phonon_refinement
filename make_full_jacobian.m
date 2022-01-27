function jacobian = make_full_jacobian(SYMS);
% Assembles the "full" Jacobian from the constituent AUX data.
%
% 

debug = 0;

n_sym = length(SYMS);
jac_nnz = SYMS{1}.VARS.jac_nnz;
n_cen = size(SYMS{1}.AUX.freevars, 1) - 1;


nQs = [];
nEs = [];
jac_rows = zeros(sum(jac_nnz), 1); % pre-declare index arrays; faster, and we can confirm total indexing at the end
jac_cols = zeros(sum(jac_nnz), 1);
jac_vals = zeros(sum(jac_nnz), 1);
jac_mask = zeros(sum(jac_nnz), 1);
prev_rows = 0;
prev_indices = 0;
for i_sym = 1:n_sym
    if debug; disp(''); disp(["Building jacobian using data from SYM : " num2str(i_sym)]); end;

    SYM = SYMS{i_sym};
    nQ = SYM.AUX.Nq;
    nE = length(SYM.AUX.eng);
    freevars = SYM.AUX.freevars;
    mask = SYM.AUX.mask;

    % append number of energy/Q points to our running total
    nQs = [nQs nQ];
    nEs = [nEs nE];

    % helpers for indexing
    nr = size(freevars, 1);
    irows = [1:nr]';
    inds_q = ones(nr, 1);
    inds_cens = ones(nr-1, 1);
    
    % loop over Q-points
    inds_jac = SYM.AUX.jac_inds;
    vars_mask = SYM.AUX.vars_mask;
    for iq = 1:nQ
        if debug; disp(["    Building jacobian using data from Q-point : " num2str(iq)]); end;

        % indices of AUX jacobian data for this Q
        aux_icens = inds_jac(find(freevars(:, 1, 1)), 1, 1)(:);
        aux_iwids = inds_jac(find(freevars(:, 1, 2)), 1, 2)(:);
        aux_ihts  = inds_jac(find(freevars(:, 1+iq, 1)), 1+iq, 1)(:);
        aux_ires  = inds_jac(find(freevars(:, 1+iq, 2)), 1+iq, 2)(:);
        aux_var = [aux_icens; aux_iwids; aux_ihts; aux_ires]';      % columns / variables
        aux_eng = SYM.AUX.eng_inds(find(mask(:, iq)), iq)(:);       % rows / energies
        aux_cols = repmat(aux_var, length(aux_eng), 1)(:);
        aux_rows = repmat(aux_eng, length(aux_var), 1)(:);
        aux_inds = sub2ind(size(SYM.AUX.auxjac), aux_rows, aux_cols);
%        if debug; disp(["length(aux_inds) :", num2str(length(aux_inds))]); end;


        % indices for full jacobian
        rows_active = prev_rows + (iq-1)*nE + find(mask(:,iq));
        vars_mask_cens = vars_mask(find(freevars(:, 1, 1)), 1, 1)(:);
        vars_mask_wids = vars_mask(find(freevars(:, 1, 2)), 1, 2)(:);
        vars_mask_hts = vars_mask(find(freevars(:, 1+iq, 1)), 1+iq, 1)(:);
        vars_mask_res = vars_mask(find(freevars(:, 1+iq, 2)), 1+iq, 2)(:);
        cols_active = [vars_mask_cens; vars_mask_wids; vars_mask_hts; vars_mask_res];
        assert(length(rows_active) == length(aux_eng), "The number of non-zero energy points must be the same in AUX and VARS")
        assert(length(cols_active) == length(aux_var), "The number of parameters being fit must be the same in AUX and VARS")


        % permute so we have every row/col combination
        cols_full = repmat(cols_active(:)', length(rows_active), 1);
        rows_full = repmat(rows_active(:), length(cols_active), 1)(:);


        % update index arrays
        jaux_vals =  full(SYM.AUX.auxjac(aux_inds));
        assert(length(jaux_vals) == length(rows_full))

        n_added = length(jaux_vals);
        ind_update = prev_indices + [1:n_added];
        jac_rows(ind_update) = rows_full;
        jac_cols(ind_update) = cols_full;
        jac_vals(ind_update) = jaux_vals;
        jac_mask(ind_update) = 1;
        prev_indices += n_added;

        if debug
            if 0 %sum(jaux_vals == 0) > 0
                i_sym
                iq
                length(jac_rows)
                length(find(jaux_vals == 0))
                surprise_zeros = (jaux_vals == 0);
%                [jac_rows(surprise_zeros) jac_cols(surprise_zeros)]
                [r, c] = ind2sub(size(SYM.AUX.auxjac), aux_inds(find(surprise_zeros)));
                [r,c]
                pause()
%                jaux_vals
            end
        end

    end  % ind Q loop

    if debug
        assert(length(find(jac_vals) == jac_nnz(i_sym)))
    end


    added_rows = nQ * nE;
    prev_rows = prev_rows + added_rows;
end % end sym loop

assert(length(jac_vals) == sum(jac_nnz), ["There's an indexing error: too many values, length " num2str(length(jac_vals)) " (expected " num2str(sum(jac_nnz)) ")."])
%assert(length(find(jac_vals)) == length(jac_vals), ["There's an indexing error: not enough values, only " num2str(length(find(jac_vals))) " (expected " num2str(sum(jac_nnz)) ")."])

n_jrows = sum(nQs .* nEs);
n_jcols = 2*n_cen +2*sum((n_cen + 1)*nQs);
jacobian = sparse(jac_rows, jac_cols, jac_vals, n_jrows, n_jcols);

if debug; 
    disp(["... Jacobian is " num2str(n_jrows) "x" num2str(n_jcols) ", with " num2str(nnz(jacobian)) " non-zero."]);

    fit_cols = sum(jacobian, 1);
    disp(["# of fittable variables in Jacobian : ", num2str(length(find(fit_cols)))]);
end;