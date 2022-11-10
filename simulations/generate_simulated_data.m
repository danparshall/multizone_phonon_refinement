function SYMS = generate_simulated_data()
addpath('../0_subs/');
addpath('../')


debug = 1;
run_fit = 1;
show_plots = 1;


%% manual, for hypothetical cubic crystal (simulation only handles cubic)
XTAL.latt = 4.287;
XTAL.cens = [6, 8, 11, 15, 20, 22, 25, 30];
%XTAL.cens = [4, 7, 10];

%XTAL.cens = [7.92   12.02  17  20.96   26.49  30  33.45  36  40.09  45  59.83  62.06] ;
XTAL.cens = [8, 12, 17, 21, 26.5, 30, 33.5, 36, 40, 45, 60, 62];
%XTAL.cens = [4, 4.5, 5, 5.5, 6, 7];
%XTAL.cens = [15, 25, 40, 4, 4.5, 5];
XTAL.wids = 0.04 * ones(size(XTAL.cens));

junk_scale = 0.2
pred_error = 0.07;   % How close the starting values of the parameters are to the true values; lower values simulate a more accurate DFT prediction 


SYMS = {};
if 1
    %% user-defined params
    SYM.Ei=100;
    SYM.chopfreq=600;
    SYM.eng=[4 : 0.25 : SYM.Ei]';
    SYM.peak_asymmetry = 1.7;
    max_Qs = 0;

    % simulated data (don't edit)
    [SYM, sim_vars] = simulate_phonon_predictions(SYM, XTAL, max_Qs);
    [SYM,startvars,true_y] = simulate_data_from_predictions(SYM, sim_vars, junk_scale, pred_error);
    SYM.sim_vars = sim_vars;
    SYM.startvars = startvars;
    SYM.true_y = true_y;
    SYM.DAT.peak_asymmetry = SYM.peak_asymmetry;
    SYM.DAT.Ei = SYM.Ei;
    SYM.DAT.chopfreq = SYM.chopfreq;
    n_sym = length(SYMS) + 1;
    SYMS{n_sym} = SYM;
end



if 0
    %% user-defined params
    SYM.Ei=50;
    SYM.chopfreq=600;
    SYM.eng=[2: 0.5 :SYM.Ei]';
    SYM.peak_asymmetry = 1.7;
    max_Qs = 50

    % simulated data (don't edit)
    [SYM, sim_vars] = simulate_phonon_predictions(SYM, XTAL, max_Qs);
    [SYM,startvars,true_y] = simulate_data_from_predictions(SYM, sim_vars, junk_scale, pred_error);
    SYM.sim_vars = sim_vars;
    SYM.startvars = startvars;
    SYM.true_y = true_y;
    SYM.DAT.peak_asymmetry = SYM.peak_asymmetry;
    SYM.DAT.Ei = SYM.Ei;
    SYM.DAT.chopfreq = SYM.chopfreq;
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


if 0
    for ind = 1:SYMS{1}.AUX.Nq
        this_dat = SYMS{1}.DAT.y_dat(SYMS{1}.AUX.mask(:,ind));
        this_err = SYMS{1}.DAT.e_dat(SYMS{1}.AUX.mask(:,ind));
        avg_delta = mean(abs(diff(this_dat, 1, 1)));
        avg_err = mean(this_err);
        disp(mean(avg_delta/avg_err))

    end
end

% fit data
if run_fit
    SYMS = refine_multizones(SYMS);

    % === plotting ===
    if show_plots
        for ind_sym = 1:length(SYMS)
            disp(['Plotting SYM ', num2str(ind_sym)]);
            SYM = SYMS{ind_sym};
            AUX = SYM.AUX;

            disp(['  ind_sym : ' num2str(ind_sym)])
            if debug
                n_cols = 12;
                disp(['Allvars size : ' num2str(size(SYMS{1}.VARS.allvars))]);
                disp('startvars :')
                SYM.startvars(:, 1:n_cols)
%                disp('Cens :')
%                cen = SYMS{1}.VARS.allvars(1:SYMS{1}.VARS.Nph)
                disp('auxvars (page1) :')
                disp(AUX.auxvars(:,1:n_cols,1))
                disp('---------------')
                disp('')

                disp('uncertainty (page1) :')
                if isfield(AUX, "uncertainty")
                    disp(AUX.uncertainty(:,:,1))
                else
                    disp("No uncertainty calculation done; edit 'refine_multizones'")
                end
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
