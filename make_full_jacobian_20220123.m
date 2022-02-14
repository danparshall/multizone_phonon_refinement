function jacobian = make_full_jacobian(SYMS);
% Assembles the "full" Jacobian from the constituent AUX data.
%
% 

debug = 1;

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
    inds_jac = SYMS{i_sym}.AUX.jac_inds;
    for iq = 1:nQ
        if debug; disp(["    Building jacobian using data from Q-point : " num2str(iq)]); end;


        % indices of jacobian data for this Q
        hts = sub2ind(size(freevars), irows, (1+iq)*inds_q, 1*inds_q);        % heights and constant BG
        res = sub2ind(size(freevars), irows, (1+iq)*inds_q, 2*inds_q);        % resolution and linear BG
        rows_aux = (iq-1)*nE + find(mask(:,iq));
        cols_aux = [ sub2ind(size(freevars), [1:nr-1]', inds_cens, inds_cens);      % centers
                    sub2ind(size(freevars), [1:nr-1]', inds_cens, 2*inds_cens);     % widths
                    hts(find(freevars(:, 1+iq, 1)));           % heights and constant BG
                    res(find(freevars(:, 1+iq, 2)));           % resolution and linear BG
                ]
        cols_jaux = repmat(cols_aux(:), length(rows_aux), 1);
        rows_jaux = repmat(rows_aux(:)', length(cols_aux), 1);
        rows_jaux = rows_jaux(:);
        inds_jaux = sub2ind(size(SYM.AUX.auxjac), rows_jaux, cols_jaux);        % get values from aux jacobian


        % indices for full jacobian
        rows_active = prev_rows + (iq-1)*nE + find(mask(:,iq));
        vars_mask_hts = SYM.AUX.vars_mask(:, 1+iq, 1);
        vars_mask_res = SYM.AUX.vars_mask(:, 1+iq, 2);
        cols_active = [[1:2*n_cen]'; vars_mask_hts(find(vars_mask_hts)); vars_mask_res(find(vars_mask_res))];
        assert(length(rows_active) == length(rows_aux), "The number of non-zero energy points must be the same in AUX and VARS")
        assert(length(cols_active) == length(cols_aux), "The number of parameters being fit must be the same in AUX and VARS")


        % permute so we have every row/col combination
        cols_add = repmat(cols_active(:)', length(rows_active), 1);
        rows_add = repmat(rows_active(:), length(cols_active), 1)(:);


        % update index arrays
        jaux_vals =  full(SYM.AUX.auxjac(inds_jaux));
        assert(length(jaux_vals) == length(rows_add))

        n_added = length(jaux_vals);
        ind_update = prev_indices + [1:n_added];
        jac_rows(ind_update) = rows_add;
        jac_cols(ind_update) = cols_add;
        jac_vals(ind_update) = jaux_vals;
        jac_mask(ind_update) = 1;
        prev_indices += n_added;

        if debug
            if sum(jac_vals(find(jac_mask))) > 0
                i_sym
                iq
                length(jac_rows)
                length(find(jac_vals == 0))
                surprise_zeros = jac_vals(find(jac_mask)) == 0;
%                [jac_rows(surprise_zeros) jac_cols(surprise_zeros)]
                [r, c] = ind2sub(size(SYM.AUX.freevars), inds_jaux)
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

assert(length(jac_vals) == sum(SYMS{1}.VARS.jac_nnz), ["There's an indexing error (too many values, length " num2str(length(jac_vals)) ")."])
assert(length(find(jac_vals)) == length(jac_vals), ["There's an indexing error (not enough values, only " num2str(length(find(jac_vals))) ")."])

n_jrows = sum(nQs .* nEs);
n_jcols = 2*n_cen +2*sum((n_cen + 1)*nQs);
jacobian = sparse(jac_rows, jac_cols, jac_vals, n_jrows, n_jcols);
