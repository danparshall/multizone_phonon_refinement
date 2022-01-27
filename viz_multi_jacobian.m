function viz_multi_jacobian()

n_atoms = 5;
n_cen = 3*n_atoms;


nQs = [ 3,  2,  4]
nEs = [20, 35, 30]
assert(length(nQs) == length(nEs), 'Need to define nQ and nE for each SYM')


n_sym = length(nQs)



%%%% ==== First loop, creating simulated data in each SYM ===
cols_prev = 2*n_cen;
jac_nnz = [];
nDatas = [];
for i_sym = 1:n_sym
    disp('')
    disp(["Generating data for SYM : ", num2str(i_sym)])
    nQ = nQs(i_sym);
    nE = nEs(i_sym);

    % make simulated mask & auxjacobian
    mask = ones(nE, nQ);    % the limiting (and impossible) case in which every datapoint contributes to fitting
    auxjac = sim_auxjac(n_cen, nQ, nE);
    if 1 %i_sym == 1
        spy(auxjac)
        disp(["   NNZ : " num2str(nnz(auxjac))])
        pause
    end
    SYMS{i_sym}.AUX.mask = mask;
    SYMS{i_sym}.AUX.auxjac = auxjac;


    % make simulated freevars
    freevars = zeros(n_cen+1, nQ+1, 2);
    freevars(:, 1, :) = 1;               % centers and widths
    freevars(:, 2:end, 1) = 1;          % assuming that all heights and constant BG can be fit
    freevars(end, 2:end, 2) = 1;        % linear BG
    freevars(end, 1, :) = 0;             % placeholder vars
    freevars(1:end-1, 2:end, 2) = 0;     % resolution widths
    SYMS{i_sym}.AUX.freevars = freevars;


%    freevars([10:15], [2:3], 1) = 0
%    freevars(end-2, 1, 2) = 0

    % 
    [vars_mask, num_nz] = make_vars_mask(SYMS{i_sym}.AUX, cols_prev);

    cols_added = 2 * nQ * (n_cen+1);
    cols_prev = cols_prev + cols_added;

    SYMS{i_sym}.AUX.vars_mask = vars_mask;
    SYMS{i_sym}.AUX.Nq = nQ;
    SYMS{i_sym}.AUX.eng = [1:nE];
    SYMS{i_sym}.AUX.num_nz = num_nz;
    jac_nnz = [jac_nnz, num_nz];
    nDatas = [nDatas nQ*nE];
end   % end first sym loop
SYMS{1}.VARS.jac_nnz = jac_nnz;
SYMS{1}.VARS.nDatas = nDatas;
disp(["num non-zero (from AUXs): ", num2str(sum(jac_nnz))])




%%%% ==== Second loop, constructing VARS.allvars and full jacobian from the SYMS ====
jacobian = make_full_jacobian(SYMS);
disp(["num non-zero (direct): ", num2str(nnz(jacobian))])
spy(jacobian)



end % end main function