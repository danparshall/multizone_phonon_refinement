function SYMS = generate_simulated_data()
addpath('../0_subs/');
addpath('../')


debug = 1;
run_fit = 1;
show_plots = 1;


%% manual, for hypothetical cubic crystal (simulation only handles cubic)
XTAL.latt = 4.287;
%XTAL.cens = [4, 7, 10, 15, 25, 30];
XTAL.cens = [4, 4.5, 5, 5.5, 6, 7];
%XTAL.cens = [15, 25, 40, 4, 4.5, 5];
XTAL.wids = 0.04 * ones(size(XTAL.cens));

junk_scale = 0.05
pred_error = 0.01;   % How close the starting values of the parameters are to the true values; lower values simulate a more accurate DFT prediction 


SYMS = {};
if 1
    %% user-defined params
    SYM.Ei=30;
    SYM.chopfreq=600;
    SYM.eng=[3 : 0.1 : SYM.Ei]';
    max_Qs = 30;

    % simulated data (don't edit)
    [SYM, sim_vars] = simulate_phonon_predictions(SYM, XTAL, max_Qs);
    sim_vars
    [SYM,startvars,true_y] = simulate_data_from_predictions(SYM, sim_vars, junk_scale, pred_error);
    SYM.sim_vars = sim_vars;
    SYM.startvars = startvars;
    SYM.true_y = true_y;
    n_sym = length(SYMS) + 1;
    SYMS{n_sym} = SYM;
end



if 1
    %% user-defined params
    SYM.Ei=50;
    SYM.chopfreq=600;
    SYM.eng=[2: 0.25 :SYM.Ei]';
    max_Qs = 50

    % simulated data (don't edit)
    [SYM, sim_vars] = simulate_phonon_predictions(SYM, XTAL, max_Qs);
    [SYM,startvars,true_y] = simulate_data_from_predictions(SYM, sim_vars, junk_scale, pred_error);
    SYM.sim_vars = sim_vars;
    SYM.startvars = startvars;
    SYM.true_y = true_y;
    n_sym = length(SYMS) + 1;
    SYMS{n_sym} = SYM;
end


% make sure starting centers agree for all conditions
for i_sym = 1:length(SYMS)
    if i_sym == 1
        startvars = SYMS{i_sym}.startvars;
    else
        SYMS{i_sym}.startvars(:,1) = startvars(:,1);
        SYMS{i_sym}.startvars(:,2) = startvars(:,2);
    end
end



%% generate AUX and VARS structures prior to fitting
SYMS = generate_AUX(SYMS);
SYMS = assemble_VARS(SYMS);

VARS = SYMS{1}.VARS;
disp([' length(y_obs) : ', num2str(length(VARS.y_obs)) ';  sum(isnan(y_obs)) : ', num2str(sum(isnan(VARS.y_obs)))]);
disp([' length(w_obs) : ', num2str(length(VARS.w_obs)) ';  sum(isnan(w_obs)) : ', num2str(sum(isnan(VARS.w_obs)))]);



% fit data
if run_fit
    SYMS = refine_multizones(SYMS);

    % === plotting ===
    if show_plots
        disp('Plotting...')
        for ind_sym = 1:length(SYMS)
            SYM = SYMS{ind_sym};
            AUX = SYM.AUX;

            disp(['  ind_sym : ' num2str(ind_sym)])
            if debug
                disp(['Allvars size : ' num2str(size(SYMS{1}.VARS.allvars))]);
                disp('startvars :')
                SYM.startvars
%                disp('Cens :')
%                cen = SYMS{1}.VARS.allvars(1:SYMS{1}.VARS.Nph)
                disp('auxvars (page1) :')
                disp(AUX.auxvars(:,:,1))
                disp('---------------')
                disp('')

                disp('uncertainty (page1) :')
                disp(AUX.uncertainty(:,:,1))
                disp('---------------')
                disp('')
            end

            % plot data, fit, and possibly true underlying function
            if 0
                plot_SYM_and_true(SYM);                 % plot just the fits
            else
                plot_SYM_and_true(SYM, SYM.true_y);     % plot the true underlying data as well
            end
        end
        hold off
    end

end  % end fit


if show_plots
    disp('Plotting jacobian structure');
    [func_out, jac_full, SYMS] = calc_full_model(SYMS);
    clf(gcf);
    spy(jac_full);
    fflush(stdout);
    pause();
end

end % endfunction