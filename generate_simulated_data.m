function SYMS = generate_simulated_data()


max_Qs = 20;

SYMS = {};

SYM.Ei=25;
SYM.chopfreq=300;
SYM.eng=[2:.25:SYM.Ei]';
[SYM, sim_vars] = simulate_phonon_predictions(SYM, max_Qs);
[SYM,startvars,true_y] = simulate_data_from_predictions(SYM, sim_vars, junk_scale=0.2);
SYM.sim_vars = sim_vars;
SYM.startvars = startvars;
SYM.true_y = true_y;
SYMS{1} = SYM;




if 1
    SYM.Ei=60;
    SYM.chopfreq=300;
    SYM.eng=[2 : 0.5 : SYM.Ei]';
    [SYM, sim_vars] = simulate_phonon_predictions(SYM, max_Qs);
    [SYM,startvars,true_y] = simulate_data_from_predictions(SYM, sim_vars, junk_scale=0.2);
    SYM.sim_vars = sim_vars;
    SYM.startvars = startvars;
    SYM.startvars(:,1) = SYMS{1}.startvars(:,1);    % make sure starting centers agree for both conditions
    SYM.startvars(:,2) = SYMS{1}.startvars(:,2);
    SYM.true_y = true_y;
    SYMS{2} = SYM;
end



%% generate AUX and VARS structures prior to fitting
SYMS = generate_AUX(SYMS);
SYMS = make_VARS(SYMS);



% fit data
SYMS = refine_phonons_multizones(SYMS);

fitcens=SYMS{1}.VARS.allvars(:,1);


end % endfunction