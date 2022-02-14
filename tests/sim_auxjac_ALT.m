function [auxjac, simple_jac] = sim_auxjac(AUX);



iter_vars = AUX.freevars;
iter_vars(:, 2:end, :) = 0;


simple_jac = [];
for iq = 1:AUX.Nq
    iter_vars(:, 1+iq, :) = AUX.freevars(:, 1+iq, :);
    simple_jac = [simple_jac; AUX.mask(:, iq) * iter_vars(:)'];
    iter_vars(:, 1+iq, :) = 0;
end
simple_jac = sparse(simple_jac);
size(simple_jac)
disp(['NNZ simple meth : ', num2str(length(find(simple_jac)))])
%spy(simple_jac)
%pause()


inds_jac = AUX.jac_inds;
auxjac_rows = [];
auxjac_cols = [];
for iq = 1:AUX.Nq
    active_cols = [ inds_jac(:, 1, 1)(:);     % centers
                    inds_jac(:, 1+iq, 1)(:);    % heights + const
                    inds_jac(:, 1, 2)(:);       % widths
    %            inds_jac(1:end-1, 1+iq, 2)(:);  % reswidth
                    inds_jac(:, 1+iq, 2)(:);  % linear BG
                    ];
    active_cols = active_cols(find(active_cols))(:)';
    active_rows = AUX.eng_inds(:, iq);
    active_rows = active_rows(find(active_rows))(:);
    jrows = repmat(active_rows, length(active_cols), 1);        % column vector with [eng; eng; eng] (i.e, iterates over energies n_col times)
    jcols = repmat(active_cols, length(active_rows), 1)(:);     % column vector with [c1; c1; c1; c2; c2; c2] (ie. has each column n_eng times)
    [length(jrows), length(jcols)]
    auxjac_rows = [auxjac_rows; jrows];
    auxjac_cols = [auxjac_cols; jcols];
end
disp(['NNZ index meth : ', num2str(length(find(auxjac_rows)))])
jvals = ones(size(auxjac_rows));
n_jrows = length(AUX.eng)*AUX.Nq
n_jcols = numel(AUX.freevars)
auxjac = sparse(auxjac_rows, auxjac_cols, jvals, n_jrows, n_jcols);
%unique(auxjac_rows)
%unique(auxjac_cols)
%auxjac([1:30])
spy(auxjac)
assert(size(simple_jac) == size(auxjac))
assert(find(simple_jac(:)) == find(auxjac(:)))
%simple_jac ~= auxjac


%assert(simple_jac == auxjac)
