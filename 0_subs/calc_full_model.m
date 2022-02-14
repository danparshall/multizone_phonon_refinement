function [func_out, jac_out, SYMS, varargout] = calc_full_model(SYMS,varsin,varargin)
%	Calculates the function and Jacobian (ensemble of partial derivatives) for
%	all Q in the SYMS set.
%
%	Size of model for a single Q is just nEng.  Size of Jacobian for a single Q
%	is nEng x (3*nPhonon (cens,wids,hts) + BGconst + BGslope)
%
%	auxjac is the jacobian for all Q-points of a given SYM instance.
%		Size is nEng * nQ rows, cols is number of elements in AUX.auxvars
%
%	jac_out is the various auxjac combined.  Ditto func_out.


debug = 0;
show_plots = 0;

if nargin > 1
	SYMS = update_AUX(SYMS,varsin);
end

Nph = SYMS{1}.VARS.Nph;

func_out = [];
for ind_sym=1:length(SYMS)

    if debug
        disp('')
        disp(['Calculating for SYM : ', num2str(ind_sym)])
        disp(['VARS.jac_nnz : ', num2str(SYMS{1}.VARS.jac_nnz)])
    end
	clear AUX;
	clear DAT;
	AUX = SYMS{ind_sym}.AUX;
	DAT = SYMS{ind_sym}.DAT;
    Nph = AUX.Nph;
    Nq = AUX.Nq;
    Ne = length(AUX.eng);

    y_calc = zeros(size(DAT.y_dat));


    aux_nnz = AUX.aux_nnz;              % calculated in make_vars_mask
    jac_vals = zeros(aux_nnz, 1);
    jac_rows = zeros(aux_nnz, 1);
    jac_cols = zeros(aux_nnz, 1);
    n_jrows = Nq * length(DAT.eng);
    n_jcols = numel(AUX.auxvars);

    % helpers for indexing
    nr = size(AUX.freevars, 1);
    curr_sparserow = 0;

    jac_inds = AUX.jac_inds;
    for iq = 1:Nq
        if debug; disp(['iq : ' num2str(iq)]); 
            AUX.auxvars(:, [1, iq+1])
        end;

        [model,jacobian] = calc_singleQ(AUX, iq);        % internally, we only evaluate at the valid energies
        valid_E = find(AUX.mask(:, iq));                 % so we will only update the model at those energies
        y_calc(valid_E, iq) = model;

        if debug; 
            disp(['jacobian shape : ', num2str(size(jacobian))]); 
            disp(['  Q NNZ : ' num2str(length(find(jacobian))) ]);

            if show_plots

                fit_eng = AUX.eng(valid_E);
                cens = AUX.auxvars(1:end-1, 1);
                hts = AUX.auxvars(1:end-1, iq+1);
                hold off;
                plot(fit_eng, model, 'r-','linewidth',1);     % best-fit line
                hold on; 
                plot(cens, hts, 'g*','linewidth',1);    % fitted centers
                
                axis([AUX.eng(1) AUX.eng(end)]);
                xlabel('Energy (meV)');
                ylabel('Intensity (arb. units)');
                legend('Calculated', 'Peaks')


                plot([0 AUX.eng(end)],[0 0],'k--');     % x-axis

        %        [fit_eng(:), model(:)]
                pause();
            end
        end;

        % most of the effort here is in constructing AUX.auxjac
        % jacobian has number of rows equal to length(valid_E); here we pick out those columns which correspond to variables being refined
        i_cens = find(AUX.freevars(:, 1, 1));
        i_wids = find(AUX.freevars(:, 1, 2)) + nr;
        i_hts = find(AUX.freevars(:, 1+iq, 1)) + 2*nr;
        i_res = find(AUX.freevars(:, 1+iq, 2)) + 3*nr;
        jvals = [reshape(jacobian(:, i_cens), [], 1);        % centers
                reshape(jacobian(:, i_hts), [], 1);          % heights + constant BG
                reshape(jacobian(:, i_wids), [], 1);         % phonon widths
                reshape(jacobian(:, i_res), [], 1);          % linear (and reswids, though not fittable)
             ];

        % indices of the rows/columns from this Q - note that these must correspond correctly to the values array, above
        active_cols = [ reshape(jac_inds(:, 1, 1), [], 1);      % centers
                        reshape(jac_inds(:, 1+iq, 1), [], 1);   % heights + const
                        reshape(jac_inds(:, 1, 2), [], 1);      % widths
                        reshape(jac_inds(:, 1+iq, 2), [], 1);   % linear BG
                        ];
        active_cols = active_cols(find(active_cols))';           % row vec containing indices of fitted columns
        active_rows = AUX.eng_inds(:, iq);
        active_rows = active_rows(find(active_rows));
        jrows = repmat(active_rows, length(active_cols), 1);     % column vector with [eng; eng; eng] (i.e, iterates over energies n_col times)
        jcols = repmat(active_cols, length(active_rows), 1);     % column vector with [c1; c1; c1; c2; c2; c2] (ie. has each column n_eng times)
        jcols = jcols(:);
        
        % update the larger array with results from this Q
        num_new_rows = length(jvals);               % new **non-zero** rows
        assert(length(jrows) == num_new_rows);
        assert(length(jcols) == num_new_rows);
        jac_rows(curr_sparserow + [1:num_new_rows]) = jrows;
        jac_cols(curr_sparserow + [1:num_new_rows]) = jcols;
        jac_vals(curr_sparserow + [1:num_new_rows]) = jvals;
        curr_sparserow = curr_sparserow + num_new_rows;

    end % end Q loop
    assert(length(jac_vals) == aux_nnz, "Indexing error, number of non-zeros is wrong")

    % attach auxjac
    SYMS{ind_sym}.AUX.auxjac = sparse(jac_rows, jac_cols, jac_vals, n_jrows, n_jcols);
    if debug && show_plots; spy(SYMS{ind_sym}.AUX.auxjac); end;

    func_out = [func_out; y_calc(:)];   % append the calculated function as a column vector; only have to do this once per SYM, so not too much overhead
end   % end SYMS loop

jac_out = make_full_jacobian(SYMS);
end     % end function