function [func_out, jac_out] = calc_full_model(SYMS,varsin,varargin)
%	Calculates the function and Jacobian (ensemble of partial derivatives) for
%	all Q in the SYMS set.
%
%	Size of model for a single Q is just nEng.  Size of Jacobian for a single Q
%	is nEng x (3*nPhonon (cens,wids,hts) + BGconst + BGslope)
%
%	jacaux is the jacobian for all Q-points of a given SYM instance.
%		Size is nEng * nQ rows, cols is number of elements in AUX.auxvars
%
%	jacout is the various jacaux combined.  Ditto funcout.


debug = 0;

if nargin > 1
	SYMS=update_AUX(SYMS,varsin);
end

Nph=SYMS{1}.VARS.Nph;

func_out = [];
for ind_sym=1:length(SYMS)

    if debug
        disp('')
        disp(['Calculating for SYM : ', num2str(ind_sym)])
    end
	clear AUX;
	clear DAT;
	AUX=SYMS{ind_sym}.AUX;
	DAT=SYMS{ind_sym}.DAT;
    Nph = AUX.Nph;
    Nq = AUX.Nq;
    Ne = length(AUX.eng);

    y_calc = zeros(size(DAT.ydat));


    jac_nnz = AUX.jac_nnz;              % calculated in make_vars_mask
%    jac_nnz = ((3*Nph) + 2) * length(find(AUX.mask));   % simple calc, assumes everything is fit except the reswidth
    jac_vals = zeros(jac_nnz, 1);
    jac_rows = zeros(jac_nnz, 1);
    jac_cols = zeros(jac_nnz, 1);
    n_jrows = Nq * length(DAT.eng);
    n_jcols = numel(AUX.auxvars);

    % helpers for indexing
    nr = size(AUX.freevars, 1);
    pagesize = numel(AUX.freevars(:,:,1));
    curr_sparserow = 0;

    inds_jac = AUX.jac_inds;
    for iq = 1:Nq
        % NOTE: if we want to include fitting resolution width, everything inside this loop should be okay; only need to implement the jacobian in "calc_singleQ"

        [model,jacobian] = calc_singleQ(AUX, iq);        % internally, we only evaluate at the valid energies
        valid_E = find(AUX.mask(:, iq));                 % so we will only update the model at those energies
        y_calc(valid_E, iq) = model;

        % jacobian has number of rows equal to length(valid_E); here we pick out those columns which correspond to variables being refined
        i_cens = find(AUX.freevars(:, 1, 1));
        i_wids = find(AUX.freevars(:, 1, 2)) + nr;
        i_hts = find(AUX.freevars(:, 1+iq, 1)) + 2*nr;
        i_res = find(AUX.freevars(:, 1+iq, 2)) + 3*nr;
        jvals = [reshape(jacobian(:, i_cens), [], 1);        % centers
                reshape(jacobian(:, i_hts), [], 1);          % heights + constant BG
                reshape(jacobian(:, i_wids), [], 1);         % phonon widths
                reshape(jacobian(:, i_res), [], 1);          % linear (and reswids, if anyone ever implements)
             ];

        % indices of the rows/columns from this Q - note that these must correspond correctly to the values array, above
        active_cols = [ reshape(inds_jac(:, 1, 1), [], 1);     % centers
                        reshape(inds_jac(:, 1+iq, 1), [], 1);    % heights + const
                        reshape(inds_jac(:, 1, 2), [], 1);       % widths
        %            reshape(inds_jac(1:end-1, 1+iq, 2), [], 1);  % reswidths
                        reshape(inds_jac(:, 1+iq, 2), [], 1);  % linear BG
                        ];
        active_cols = active_cols(find(active_cols))';           % row vec containing indices of fitted columns
        active_rows = AUX.eng_inds(:, iq);
        active_rows = active_rows(find(active_rows));
        jrows = repmat(active_rows, length(active_cols), 1);     % column vector with [eng; eng; eng] (i.e, iterates over energies n_col times)
        jcols = repmat(active_cols, length(active_rows), 1);     % column vector with [c1; c1; c1; c2; c2; c2] (ie. has each column n_eng times)
        jcols = jcols(:);
        
        % update the larger array with results from this Q
%        num_new_rows = length(valid_E) * length(active_cols);
        num_new_rows = length(jvals);
        assert(length(jrows) == num_new_rows);
        assert(length(jcols) == num_new_rows);
        jac_rows(curr_sparserow + [1:num_new_rows]) = jrows;
        jac_cols(curr_sparserow + [1:num_new_rows]) = jcols;
        jac_vals(curr_sparserow + [1:num_new_rows]) = jvals;
%        row_offset += Ne;
%        col_offset += nr;
        curr_sparserow = curr_sparserow + num_new_rows;

        if 0 % iq == 1
            [jrows jcols jvals]
        end
    end % end Q loop
%    disp(['crude : ', num2str(length(crude_vals))])
    assert(length(jac_vals) == jac_nnz, "Indexing error, number of non-zeros is wrong")
%    jac_rows
    auxjac = sparse(jac_rows, jac_cols, jac_vals, n_jrows, n_jcols);
%    size(auxjac)
%    auxjac([1:30])
%    auxjac(find(auxjac)(1:30))
    SYMS{ind_sym}.AUX.auxjac = auxjac;
    spy(SYMS{ind_sym}.AUX.auxjac)
%    [simjac, simple_jac] = sim_auxjac_ALT(AUX);
%    find(simjac) == find(SYMS{ind_sym}.AUX.auxjac)

    func_out = [func_out; y_calc(:)];
end   % end SYMS loop


full_jacobian = make_full_jacobian(SYMS);
jac_out = full_jacobian(:, SYMS{1}.VARS.indfree);

end     % end function