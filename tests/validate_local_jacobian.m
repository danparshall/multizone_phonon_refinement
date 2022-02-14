function validate_local_jacobian(SYMS);


DELTA = 0.00000001;
TOLERANCE = 0.0001;

for i_sym = 1:length(SYMS)
    SYM = SYMS{i_sym};
    AUX = SYM.AUX;

    n_vars = prod(size(AUX.freevars));
    free_idx = AUX.freevars .* reshape([1:n_vars], size(AUX.freevars));


	vsize = size(AUX.freevars, 1);
    fit_cens = find(AUX.freevars(:, 1, 1));
    fit_wids = find(AUX.freevars(:, 1, 2));

    for iq = 1:AUX.Nq
        [modelout,jacout] = calc_singleQ(AUX,iq);
        jac_empirical = zeros(size(jacout));


        %% centers
        for ic = 1:length(fit_cens)
            icen = fit_cens(ic);
            TMPAUX = AUX;
            cenval = TMPAUX.auxvars(icen, 1, 1);
            
            TMPAUX.auxvars(icen, 1, 1) = cenval + DELTA;
            model_pos = calc_singleQ(TMPAUX, iq);

            TMPAUX.auxvars(icen, 1, 1) = cenval - DELTA;
            model_neg = calc_singleQ(TMPAUX, iq);
     
            jac_empirical(:, icen) = (model_pos - model_neg)/(2*DELTA);            
        end % cens


        %% widths
        for iw = 1:length(fit_wids)
            iwid = fit_wids(iw);
            TMPAUX = AUX;
            widval = TMPAUX.auxvars(iwid, 1, 2);
            
            TMPAUX.auxvars(iwid, 1, 2) = widval + DELTA;
            model_pos = calc_singleQ(TMPAUX, iq);

            TMPAUX.auxvars(iwid, 1, 2) = widval - DELTA;
            model_neg = calc_singleQ(TMPAUX, iq);

            jac_empirical(:, vsize + iwid) = (model_pos - model_neg)/(2*DELTA);     
        end  % wids


        %% heights and constant BG
        fit_hts = find(AUX.freevars(:, 1+iq, 1));
        for ih = 1:length(fit_hts)
            iht = fit_hts(ih);
            TMPAUX = AUX;
            htval = TMPAUX.auxvars(iht, 1+iq, 1);
            
            TMPAUX.auxvars(iht, 1+iq, 1) = htval + DELTA;
            model_pos = calc_singleQ(TMPAUX, iq);

            TMPAUX.auxvars(iht, 1+iq, 1) = htval - DELTA;
            model_neg = calc_singleQ(TMPAUX, iq);
  
            jac_empirical(:, 2*vsize + iht) = (model_pos - model_neg)/ (2*DELTA);     
        end  % hts


        %% linear BG and resolution
        fit_res = find(AUX.freevars(:, 1+iq, 2))
        for ir = 1:length(fit_res)
            ires = fit_res(ir);
            TMPAUX = AUX;
            resval = TMPAUX.auxvars(ires, 1+iq, 2);
            
            TMPAUX.auxvars(ires, 1+iq, 2) = resval + DELTA;
            model_pos = calc_singleQ(TMPAUX, iq);

            TMPAUX.auxvars(ires, 1+iq, 2) = resval - DELTA;
            model_neg = calc_singleQ(TMPAUX, iq);

            jac_empirical(:, 3*vsize + ires) = (model_pos - model_neg)/ (2*DELTA);     
        end  % res


        %% compare
        jac_diff = jacout - jac_empirical;
        if sum(sum(abs(jac_diff) > TOLERANCE)) > 0
            disp(['PROBLEM in SYM ' num2str(i_sym) ', Q-point : ', num2str(iq)])
            AUX.freevars
            [jacout, jac_empirical, jac_diff]
            pause
        end  % jac inspection

    end % q loop

end  % SYM loop


end % end function