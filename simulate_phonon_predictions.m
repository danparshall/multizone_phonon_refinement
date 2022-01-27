function [SYM, sim_vars] = simulate_phonon_predictions(SYM, max_Qs)


    function eng=calc_mom_to_eng(mom);
        % calculates neutron energy (meV) from momentum (inv Angstroms)
        eng = 81.82*(mom./(2*pi)).^2;
    end


    function mom=calc_eng_to_mom(eng);
        % Calculates neutron momentum (inv Angstroms) from energy (meV)
        mom = 2*pi*sqrt(eng/81.82);
    end


    function k_f_max = calc_kfmax(Q_mag, k_i)
        if Q_mag < k_i;
            k_f_max = k_i - Q_mag;
        else
            k_f_max = Q_mag - k_i;
        end
    end

%% manual, for hypothetical cubic crystal
latt = 4.287;
cens = [5, 7, 10, 15, 30, 31];
cens = cens(:);
cens = cens + 0.00005;

E_i = SYM.Ei;
N_ph = length(cens);
wids = 0.04 * ones(N_ph, 1);


k_i = calc_eng_to_mom(E_i);
Q_max = calc_eng_to_mom(E_i);
Q_latt = 2*pi / latt;
max_bzone = floor(Q_max / Q_latt);


HKL_vals = [];
Q_mags = [];
mask = [];
E_maxes = [];

for H = 0:max_bzone
    for K = 0:max_bzone
        for L = 0:max_bzone
            poss = [H K L];
            Q_mag = sqrt(sum((Q_latt * poss).^2));

            k_f_max = calc_kfmax(Q_mag, k_i);

            E_max = calc_mom_to_eng(k_f_max);
            poss_phonons = cens < (E_max + 2);

            if (Q_mag < Q_max) && (Q_mag > 0) && sum(poss_phonons) > 0
                HKL_vals = [HKL_vals; poss];
                Q_mags = [Q_mags; Q_mag];
                mask = [mask, poss_phonons];
                E_maxes = [E_maxes E_max];
            end  % end-if
        end
    end
end


rand('seed', 42);               % set random seed, for data reproducibility

% use only up to max_Qs points
if length(Q_mags) > max_Qs
    rand_idx = randperm(length(Q_mags))(1:max_Qs);
    HKL_vals = HKL_vals(rand_idx, :);
    Q_mags = Q_mags(rand_idx);
    mask = mask(:, rand_idx);
end

% simulate peak heights
heights = 2*rand(size(mask));   % average height is 1
heights = heights .* mask;
sim_vars = [cens(:), wids(:), heights];


SYM.DAT.HKL_vals = HKL_vals;
SYM.DAT.Q_mag = Q_mags;
SYM.E_maxes = E_maxes;  % kinematic contstraint; not needed for real data.
end  % end-function