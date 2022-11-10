function [SYM, sim_vars] = simulate_phonon_predictions(SYM, XTAL, max_Qs, seed)


% set random seed, for data reproducibility
if ~exist('seed'); seed=42; end
rand('seed', seed);


%% helper functions
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


function symm_qs = apply_cubic_ops(q_vec)
    % apply cubic symmetry operations
    symm_qs = [
        unique(perms(q_vec), 'rows');
        -1*unique(perms(q_vec), 'rows');
    ];
    symm_qs(symm_qs == 0) = 0;
end



%% only designed for cubic crystal
latt = XTAL.latt;
cens = XTAL.cens(:);
cens = cens + 0.00005;      % tiny offset, reduces pathological cases of vanishing derivatives when a center is exactly equal to an energy point
wids = XTAL.wids(:);


E_i = SYM.Ei;
N_ph = length(cens);


k_i = calc_eng_to_mom(E_i);
Q_max = calc_eng_to_mom(E_i);
Q_latt = 2*pi / latt;
max_bzone = floor(Q_max / Q_latt);


HKL_vals = [];
Q_mags = [];
mask = [];
E_maxes = [];
q_vec = [0, 0, 0.2];
disp(["Generating data for reduced q : [" num2str(q_vec(1)) ", " num2str(q_vec(2)) ", " num2str(q_vec(3)) "]"])
cubic_ops = apply_cubic_ops(q_vec);

for H = 0:max_bzone
    for K = 0:max_bzone
        for L = 0:max_bzone
            possibles = [H K L] + cubic_ops;
            for P = 1:size(possibles, 1)
                poss = possibles(P,:);
                if ismember(poss, HKL_vals, 'rows')
                    continue
                end
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
            end  % end perms
        end
    end
end


% use only up to max_Qs points
if (max_Qs > 0) && (length(Q_mags) > max_Qs)
    disp(["Restricting to " num2str(max_Qs) " out of " num2str(length(Q_mags)) " Q-points."])
    rand_idx = randperm(length(Q_mags));
    rand_idx = rand_idx(1:max_Qs);
    HKL_vals = HKL_vals(rand_idx, :);
    Q_mags = Q_mags(rand_idx);
    mask = mask(:, rand_idx);
end

% simulate peak heights
heights = rand(size(mask));   % average height is 1
heights = heights .* mask;
sim_vars = [cens(:), wids(:), heights];


%SYM.DAT.HKL_vals = HKL_vals;
SYM.DAT.Q_hkl = HKL_vals;
SYM.DAT.Qmags = Q_mags;
SYM.E_maxes = E_maxes;  % kinematic constraint; not needed for real data.
end  % end-function
