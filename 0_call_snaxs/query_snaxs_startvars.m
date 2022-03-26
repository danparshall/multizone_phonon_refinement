function [centers, widths, heights] = query_snaxs_startvars(PAR, Q_hkl, symm, varargin);
% Call SNAXS for a set of Qs, applying symmetry operators to handle twinning

if nargin == 2
    % simple case, no twinning
    [centers, widths, heights] = query_snaxs(PAR, Q_hkl);

elseif nargin == 3
    % when a symm object is passed, calculate phonons in the BZ produced by twinning (as well as the nominal BZ)

    for iq = 1:length()
        Q_twinned = symmops(Q_hkl(iq, :), symm);
        [cens, wids, hts] = query_snaxs(PAR, Q_twinned);

        if iq == 1
            centers = cens;
            widths = wids;
            heights = hts;
        else
            TOL_CENS = 1E-9;
            TOL_WIDS = 1E-6;
            assert(isclose(centers, cens, TOL_CENS), "Twinned phonon centers out of tolerance")
            assert(isclose(widths, wids, TOL_WIDS), "Twinned phonon widths out of tolerance")
            disp(['SIZE(hts) : ', num2str(size(hts))])
            heights = [heights, sum(hts, 2)];
        end  % 
    end
else
    errror
end




function [centers, widths, heights] = query_snaxs(PAR, Q_hkl)

    % === initialize arrays, index variables ===
    [XTAL,EXP,INFO,PLOT,DATA,VECS] = params_fetch(PAR);
    INFO.Q = Q_hkl(1, :);
    DATA = make_DATA(PAR);
    INFO.Q_npts = size(Q_hkl,1);
    [unique_tau, cellarray_qs, Q_hkl, Q_delta] = generate_tau_q_from_Q(PAR,Q_hkl);


    % === generate VECS (including structure factor) ===
    PAR = params_update(XTAL,EXP,INFO,PLOT,DATA,VECS);
    PAR = simulate_multiQ(PAR, Q_hkl);

    TOL_ENG = 1E-9;
    assert(isclose(PAR.VECS.energies(:,1), mean(PAR.VECS.energies, 2), TOL_ENG));  % sanity-check that energies truly are the same at all Q


    % === get energy and widths ===
    [centers, iWid] = unique(PAR.VECS.energies);
    widths = PAR.VECS.phWidths(iWid);	% each center has a unique intrinsic linewidth


    % === calculate height from structure factor ===
    hts_array = calc_height_multiQ(PAR, Q_hkl);
    height_nan = ~isnan(hts_array);


    heights = zeros(length(centers), INFO.Q_npts);
    for cen = 1:length(centers)
        i_cen = find(PAR.VECS.energies(:, 1) == centers(cen));
        heights(cen, :) = sum(hts_array(i_cen, :), 1);
    end
end % end query

end % end full function