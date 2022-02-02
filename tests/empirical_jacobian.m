function jac_empirical = empirical_jacobian(AUX, iq, DELTA);



vsize = size(AUX.freevars, 1);
fit_cens = find(AUX.freevars(:, 1, 1));
fit_wids = find(AUX.freevars(:, 1, 2));


[modelout,jacout]=calc_singleQ(AUX,iq);
jac_empirical = zeros(size(jacout));


for ic = 1:length(fit_cens)
    icen = fit_cens(ic);
    TMPAUX = AUX;
    cenval = TMPAUX.auxvars(icen, 1, 1);
    
    TMPAUX.auxvars(icen, 1, 1) = cenval + DELTA;
    model_pos = calc_singleQ(TMPAUX, iq);

    TMPAUX.auxvars(icen, 1, 1) = cenval - DELTA;
    model_neg = calc_singleQ(TMPAUX, iq);

%            jac_empirical(:, free_idx(icen, 1, 1)) = (model_pos - model_neg)/2;       
    jac_empirical(:, icen) = (model_pos - model_neg)/(2*DELTA);            
end % cens


for iw = 1:length(fit_wids)
    iwid = fit_wids(iw);
    TMPAUX = AUX;
    widval = TMPAUX.auxvars(iwid, 1, 2);
    
    TMPAUX.auxvars(iwid, 1, 2) = widval + DELTA;
    model_pos = calc_singleQ(TMPAUX, iq);

    TMPAUX.auxvars(iwid, 1, 2) = widval - DELTA;
    model_neg = calc_singleQ(TMPAUX, iq);

#            jac_empirical(:, free_idx(iwid, 1, 2)) = (model_pos - model_neg)/2; 
    jac_empirical(:, vsize + iwid) = (model_pos - model_neg)/(2*DELTA);     
end  % wids


#        fit_hts = free_idx(find(free_idx(:, 1+iq, 1)));
fit_hts = find(AUX.freevars(:, 1+iq, 1));
for ih = 1:length(fit_hts)
    iht = fit_hts(ih);
    TMPAUX = AUX;
    htval = TMPAUX.auxvars(iht, 1+iq, 1);
    
    TMPAUX.auxvars(iht, 1+iq, 1) = htval + DELTA;
    model_pos = calc_singleQ(TMPAUX, iq);

    TMPAUX.auxvars(iht, 1+iq, 1) = htval - DELTA;
    model_neg = calc_singleQ(TMPAUX, iq);

%            jac_empirical(:, free_idx(iht, 1+iq, 1)) = (model_pos - model_neg)/2;     
    jac_empirical(:, 2*vsize + iht) = (model_pos - model_neg)/ (2*DELTA);     
end  % hts


%        fit_res = free_idx(find(free_idx(:, 1+iq, 2)));
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

