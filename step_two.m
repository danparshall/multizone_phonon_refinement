function SYMS = step_two();


show_plots = 0;
perform_fit = 1;

addpath('0_call_snaxs/');
addpath('0_subs/');
path_PRED = './0_call_snaxs/PREDS_test.mat'

path_EXTR = '/home/dan/code/bkbo/data/condensed_EXTRACTED__HKL__H_0.500__K_0.500__L_0.500__T_300K____20140525T161403_.mat'
path_EXTR = '~/code/bkbo/data/condensed_EXTRACTED__HKL__H_0.300__K_0.000__L_0.000__T_300K____20140525T174525_.mat'
path_EXTR = '/home/dan/code/bkbo/data/condensed_EXTRACTED__HKL__H_0.000__K_0.000__L_0.000__T_300K____20140525T153052_.mat'

path_EXTR = '~/code/bkbo/data/condensed_EXTRACTED__HKL__H_0.500__K_0.000__L_0.000__T_300K____20140525T180025_.mat'


CENSCALE = 1;
SYMS = load_data_and_preds(path_EXTR, path_PRED, CENSCALE);

    SYMS{1}.Ei=100;
    SYMS{1}.chopfreq=600;

if 0
    % GAMMA (used censcale = 0.8)
    disp("Over-riding calcualted peak centers")
    SYMS{1}.startvars(:, 1) = [0, 6, 14, 18, 25,  60];

    disp("Over-riding peak intensities")
    SYMS{1}.startvars(2:end, 3:end) = 0.5;

    disp('over-riding phonon widths')
    SYMS{1}.startvars(:, 2) = 2;
elseif 1
    % Z-point 
    z_cens = [0. 10. 17. 23. 28. 34. 40. 44. 53. 59. 70]';
    SYMS{1}.startvars(:, 1) = z_cens;

    e_scaled = 10 ./ z_cens(2:end);
    Nq = size(SYMS{1}.startvars, 2) - 2;
%    e_scaled = [2 e_scaled];
    SYMS{1}.startvars(2:end, 3:end) = repmat(e_scaled(:), 1, Nq);

    disp('over-riding phonon widths')
    SYMS{1}.startvars(:, 2) = 1;
end


%% generate AUX and VARS structures prior to fitting
SYMS = generate_AUX(SYMS);
SYMS = assemble_VARS(SYMS);

VARS = SYMS{1}.VARS;
disp([' length(y_obs) : ', num2str(length(VARS.y_obs)) ';  sum(isnan(y_obs)) : ', num2str(sum(isnan(VARS.y_obs)))]);
disp([' length(w_obs) : ', num2str(length(VARS.w_obs)) ';  sum(isnan(w_obs)) : ', num2str(sum(isnan(VARS.w_obs)))]);



if perform_fit
    SYMS = refine_multizones(SYMS);
else
    disp("")
    disp(" Skipping fitting step")
end

if show_plots;
    for i_sym = 1:length(SYMS)
        SYM = SYMS{i_sym};
        plot_SYM(SYM);
    end
end;
