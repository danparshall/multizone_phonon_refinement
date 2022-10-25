function jacobian = manual_jac_check(SYMS);
% Confirm that the Jacobian is constructed correctly, by comparing each column of the calculated jacobian 
% to an "empirical jacobian", obtained by directly calcuating the finite difference of each variable.
% Each variable is checked against a tolerance, and if the derivatives don't match, the index of the 
% variable in the original AUX is shown, to aid in troubleshooting.
%
% Can also plot the spy() of the calculated and empirical jacobians, as well.  Note that the empirical
% jacobian falls to zero faster than the calculated version (I think because of floating-point effects),
% so these don't line up perfectly.

DELTA = 0.000001;  % Change of input variables; should be large enough to avoid floating-point errors; ~e-6 or so
TOLERANCE = 0.001; % allowed difference between the closed-form and finite-difference results; should be larger than DELTA
show_plots = 1;
debug = 1;


VARS = SYMS{1}.VARS;
[func_out, jac_out] = calc_full_model(SYMS,VARS.varsin);

if show_plots;
    fig_jac = figure('Name', 'JacFull');
    spy(jac_out);
    disp('Assembled jac plotted...')
    fflush(stdout);
    pause();
    fig_vars = figure('Name', 'Derivative');
end;

n_cen = VARS.Nph;
nQs = VARS.nQs;
nEs = VARS.nEs;


switch_cols = 2*n_cen + 2*(n_cen+1)*[0 nQs];  % indices at which columns transition to each SYM (first 2*n_cen are common to all AUX, but after it depends on nQ)
xfunc = [1:sum(nQs .* nEs)];

row_inds = cell();
col_inds = cell();
jac_vals = cell();
jrows = [];
jcols = [];
jvals = [];

i_cells = 0;
indfree = VARS.indfree;         %% varsin is only the variables that are fitted; indfree are their indices
for i = 1:length(indfree)
    if mod(i, 10) == 0
        disp(["Checking variable " num2str(i)])
    end
    ivar = indfree(i);          %% column number of jacobian for this variable
    varsin = VARS.varsin;
    var_value = varsin(i);
    if debug; [i, ivar, var_value]; end;

    varsin(i) = var_value + DELTA;
    func_pos = calc_full_model(SYMS, varsin);

    varsin(i) = var_value - DELTA;
    func_neg = calc_full_model(SYMS, varsin);

    func_deriv = full((func_pos(:) - func_neg(:))) / (2*DELTA);
    calc_deriv = jac_out(:, ivar);
    assert(length(calc_deriv) == length(func_deriv))
%    assert(find(func_deriv) == find(calc_deriv))           %% doesn't work, empirical derivative goes to zero too quickly
    diffs = calc_deriv - func_deriv;
    out_of_tol =  abs(diffs) > TOLERANCE;
    if sum(out_of_tol) > 0
        disp(['Problem with Full jacobian, column ', num2str(ivar)])
        if debug; [func_deriv, full(jac_out(:, ivar))]; end;

        % === plotting ===
        if show_plots
            if ivar <= n_cen;
                disp([' cen : ' num2str(ivar) ])
            elseif ivar <= 2*n_cen;
                disp([' wid : ', num2str(ivar - n_cen)]);
            else
                run = 1;
                for i_sym = 1:length(switch_cols);
                    if ivar < switch_cols(1);
                        disp('Error trying to map VARSIN to AUX Jacobian')
                    elseif (ivar > switch_cols(i_sym)) && (run == 1);
                        AUX = SYMS{i_sym}.AUX;
                        if i_sym == 1
                            pre_cols = 2*n_cen
                        else
                            pre_cols = 2*n_cen + 2*(n_cen+1)*sum(nQs(1:i_sym-1))
                        end;
                        page_size = (nQs(i_sym)+1) * (n_cen+1);
                        avar =  (ivar - pre_cols + n_cen+1);
                        if avar > page_size
                            % second page; res + linear, need to add in extra for widths
                            avar = avar + n_cen + 1;
                        end
                        disp(['avar : ', num2str(avar)]);
                        run == 0;
                    end;
                end;
                [pos_row, pos_col, pos_page] = ind2sub(size(AUX.freevars), avar);
                disp(['auxvars loc : ', num2str([pos_row, pos_col, pos_page])]);
                disp([' auxjac ind : ', num2str(AUX.jac_inds(pos_row, pos_col, pos_page))]);
            end;


            hold off;
            plot(xfunc, func_deriv, 'g-','linewidth',2);    % empirical jac
            hold on;
            plot(xfunc, full(jac_out(:,ivar)),'r:','linewidth',2);  % calculated jac
            legend('Empirical', 'Calculated')

            plot([1 xfunc(end)],[0 0],'k--');               % x-axis
            hold off;

            pause();
        end
    end



    rows_active = find(func_deriv);
    n_rows_active = length(rows_active);

    if n_rows_active > 0            %% cell2mat is too stupid to handle empty arrays
%        disp('appending...')
        i_cells = i_cells + 1;
        cols_active = repmat(ivar, length(rows_active), 1);
%        row_inds{i_cells} = rows_active';
%        col_inds{i_cells} = cols_active';
%        jac_vals{i_cells} = func_deriv(find(func_deriv));
        jrows = [jrows; rows_active(:)];
        jcols = [jcols; cols_active(:)];
        jvals = [jvals; func_deriv(find(func_deriv))(:)];
    end
end
%row_inds
%col_inds
jrows = cell2mat(row_inds);
jcols = cell2mat(col_inds);
jvals = cell2mat(jac_vals);


n_jrows = sum(VARS.nQs .* VARS.nEs);
n_jcols = 2*n_cen +2*sum((n_cen + 1)*VARS.nQs);
jacobian = sparse(jrows, jcols, jvals, n_jrows, n_jcols);
spy(jacobian);
end