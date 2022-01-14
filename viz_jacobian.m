function viz_jacobian()

n_atoms = 5;
n_cens = 3*n_atoms;
n_Qs = 4;
eng = [1:20];
n_eng = length(eng);


auxvars = zeros(n_cens+1, n_Qs+1, 2);
auxvars(:, 1, :) = 1;               % centers and widths
auxvars(end, 1, :) = 0;             % placeholder vars
auxvars(1:end-1, 2:end, 2) = 0;     % resolution widths
%auxvars
disp(['auxvars : ', num2str(size(auxvars))])

nnz = n_Qs * (3*n_cens + 2) * n_eng;
n_jrows = (n_Qs*n_eng);
n_jcols = numel(auxvars);
disp(['Final jacobian should have size ' num2str(n_jrows) 'x' num2str(n_jcols) ', and ' num2str(nnz) ' nonzero'])


jac_rows = [];
jac_cols = [];

for i_Q = 1:n_Qs
    startrow = ((i_Q-1) * n_eng) + 1;
    rows_active = [startrow : startrow + n_eng - 1]';
%    disp('')
%    disp(['num rows : ', num2str(length(rows_active))])

    % get active columns
    auxvars( :,  i_Q + 1, 1) = 1;   % heights + constant BG
    auxvars(end, i_Q + 1, 2) = 1;   % linear BG
    cols_active = find(auxvars);
    auxvars(:, i_Q+1, :) = 0;       % reset auxvars

    if 0 %i_Q == 1
        disp(cols_active)
        disp(['num cols : ', num2str(length(cols_active))])
    end
    for i_var = 1:length(cols_active)
        col_ind = cols_active(i_var);
        jac_rows = [jac_rows; rows_active];
        jac_cols = [jac_cols;  col_ind * ones(length(rows_active), 1)];
%        disp(['len(rows) : ' num2str(length(jac_rows))])
%        disp(['len(cols) : ' num2str(length(jac_cols))])
    end % end var loop
end  % end Q loop
%[jac_rows, jac_cols]

% shape should end up with (n_Qs * n_eng) rows, and (numel(auxvars)) columns.
% for each Q, we should have n_eng*(3*n_cens + 2) columns that are nonzero
vals = ones(nnz,1);
density = nnz / (n_jrows * n_jcols)

jacobian = sparse(jac_rows, jac_cols, vals, n_jrows, n_jcols);
spy(jacobian)
end % end-function