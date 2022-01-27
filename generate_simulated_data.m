function SYMS = generate_simulated_data()

addpath('0_subs/');
if system_octave
    pkg load optim
end

max_Qs = 20;
junk_scale = 0.001

SYMS = {};

    SYM.Ei=40;
    SYM.chopfreq=300;
    SYM.eng=[2 : 0.5 : SYM.Ei]';
[SYM, sim_vars] = simulate_phonon_predictions(SYM, max_Qs);
[SYM,startvars,true_y] = simulate_data_from_predictions(SYM, sim_vars, junk_scale);
SYM.sim_vars = sim_vars;
SYM.startvars = startvars;
SYM.true_y = true_y;
SYMS{1} = SYM;




if 1

SYM.Ei=30;
SYM.chopfreq=300;
SYM.eng=[2:.25:SYM.Ei]';

    [SYM, sim_vars] = simulate_phonon_predictions(SYM, max_Qs);
    [SYM,startvars,true_y] = simulate_data_from_predictions(SYM, sim_vars, junk_scale);
    SYM.sim_vars = sim_vars;
    SYM.startvars = startvars;
    SYM.startvars(:,1) = SYMS{1}.startvars(:,1);    % make sure starting centers agree for both conditions
    SYM.startvars(:,2) = SYMS{1}.startvars(:,2);
    SYM.true_y = true_y;
    SYMS{2} = SYM;
end



%% generate AUX and VARS structures prior to fitting
SYMS = generate_AUX(SYMS);
SYMS = assemble_VARS(SYMS);


% for i_sym = 1:length(SYMS)   % this is worst-case scenario, since it implies every height impacts every eng
%     AUX = SYMS{i_sym}.AUX;
%     auxjac = AUX.mask(:) * AUX.freevars(:)';
%     disp(['NNZ : ', num2str(length(find(auxjac)))])
%     SYMS{i_sym}.AUX.auxjac = auxjac;
% end  % end auxjac simulation loop


% fit data
if 1
    SYMS = refine_multizones(SYMS);
%    fitcens=SYMS{1}.VARS.allvars(:,1);
end  % end fit

end % endfunction