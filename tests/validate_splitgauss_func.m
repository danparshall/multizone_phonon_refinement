function validate_splitgauss_func()

eng = [-2: 0.1 : 2];
cen = 0;
ht = 1;

DELTA = 0.00000001;
TOLERANCE = 0.0001;


if 0
    hwhm_1 = 1;
    hwhm_2 = 2;
    [model, jac_out] = calc_splitgauss_JAC_fast(eng, cen, ht, hwhm_1, hwhm_2);
    model_pos = calc_splitgauss_JAC_fast(eng, cen, ht, hwhm_1 + DELTA, hwhm_2 + DELTA);
    model_neg = calc_splitgauss_JAC_fast(eng, cen, ht, hwhm_1 - DELTA, hwhm_2 - DELTA);
else
    wid_ph = 0.4;
    wid_res = 0.8;
    asymm = 1.7;
    [model, jac_out] = calc_splitgauss_full(eng, cen, ht, wid_ph, wid_res, asymm);
    model_pos = calc_splitgauss_full(eng, cen, ht, wid_ph + DELTA, wid_res, asymm);
    model_neg = calc_splitgauss_full(eng, cen, ht, wid_ph - DELTA, wid_res, asymm);

end   % end if

jac_emp = (model_pos - model_neg) / (2*DELTA);

jac_diff = jac_out(:, 3) - jac_emp;


[jac_out(:,3), jac_emp, jac_diff, jac_out(:,3)./jac_emp]
if sum(abs(jac_diff) > TOLERANCE) > 0
    disp("ERROR in splitgauss jacobian")
end