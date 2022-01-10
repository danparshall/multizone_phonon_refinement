function van=get_chop_DAN(ei,omega,e_trans)
global x0 xa x1 x2 wa_mm ha_mm wa ha pslit dslat radius rho tjit mod_type s thetam mod_type imod sx sy sz isam gam ia ix idet dd tbin titledata
convert = 2.354820d0;

[v_van, v_van_m, v_van_ch, v_van_jit, v_van_ya, v_van_x, v_van_y, v_van_xy, v_van_dd]=van_var(ei,omega);
%do the conversions for the resolutions

van = real( convert .* 8.747832d-4 .* sqrt(ei.^3) .* ( sqrt(v_van + (tbin.^2./12.0d0)) .* 1.0d6 ) ./ x2 );

for i=1:length(e_trans)
	etrans=e_trans(i);
    [v_van, v_van_m, v_van_ch, v_van_jit, v_van_ya, v_van_x, v_van_y, v_van_xy, v_van_dd]=van_var(ei,omega,etrans);
    van(i) = real( convert .* 8.747832d-4 .* sqrt((ei-etrans).^3) .* ( sqrt(v_van + (tbin.^2./12.0d0)) .* 1.0d6 ) ./ x2 );
end
