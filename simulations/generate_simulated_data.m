function SYMS = generate_simulated_data()
addpath('../0_subs/');
addpath('../')


debug = 1;
run_fit = 1;
show_plots = 1;


%% manual, for hypothetical cubic crystal (simulation only handles cubic)
XTAL.latt = 4.287;
XTAL.cens = [5, 6, 7, 15, 25, 30];
XTAL.wids = 0.04 * ones(size(XTAL.cens));

junk_scale = 0.05
pred_error = 0.07;   % How close the starting values of the parameters are to the true values; lower values simulate a more accurate DFT prediction 


SYMS = {};
if 1
    %% user-defined params
    SYM.Ei=60; 
    SYM.chopfreq=300;
    SYM.eng=[2 : 0.5 : SYM.Ei]';
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
    SYM.Ei=40;
    SYM.chopfreq=600;
    SYM.eng=[2: 0.25 :SYM.Ei]';
    max_Qs = 30

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

% fit data
if run_fit
    SYMS = refine_multizones(SYMS);
%    fitcens=SYMS{1}.VARS.allvars(:,1);


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
                disp('auxvars :')
                disp(AUX.auxvars)
                disp('-------------\n')
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




end % endfunction