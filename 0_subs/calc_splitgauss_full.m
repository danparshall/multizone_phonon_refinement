function [splitgauss,jacobian,varargout]=calc_splitgauss_full(eng, cen, ht, wid_ph, wid_res, ratio_res);
% [splitgauss,jacobian]=calc_splitgauss_full(eng, cen, ht, wid_ph, wid_res, ratio_res);
% 	quickest (though confusing) method to calculate splitgauss.
% 	calculates jacobian when nargout==2
%	splitgauss size is [ len(eng) ]
%	jacobian size is [ len(eng) x 3 ]

eng = eng(:);


%%%% WIDTHS (assuming convolution of Gaussian for phonon and Gaussion for resolution; overall width is added in quadrature)
FWHM = sqrt(wid_ph^2 + wid_res^2);
con_lo = ratio_res/(ratio_res + 1);
con_hi = 1/(ratio_res + 1);
hwhm_1 = FWHM * con_lo;
hwhm_2 = FWHM * con_hi;

%splitgauss
splitgauss= [exp( -log(2)* ((eng(eng<=cen)-cen)/hwhm_1).^2);
		 	 exp( -log(2)* ((eng(eng>cen)-cen)/hwhm_2).^2)];

% jacobian
if nargout==2
	jac_cen= ht*2*log(2) * splitgauss .* [(eng(eng<=cen)-cen)/hwhm_1^2; (eng(eng>cen)-cen)/hwhm_2^2];
	jac_ht = splitgauss;
%	jac_wid= ht*2*log(2) * splitgauss .* [(eng(eng<=cen)-cen).^2./hwhm_1.^3; (eng(eng>cen)-cen).^2./hwhm_2.^3];
	jac_wid = ht*2*log(2) * wid_ph * splitgauss .* [(con_lo^2)*(eng(eng<=cen)-cen).^2/hwhm_1^4; (con_hi^2)*(eng(eng>cen)-cen).^2/hwhm_2^4];

%%%% if resolution width gets refined, need to uncomment this
%	jac_res = ht*2*log(2) * wid_res * splitgauss .* [(con_lo^2)*(eng(eng<=cen)-cen).^2/hwhm_1^4; (con_hi^2)*(eng(eng>cen)-cen).^2/hwhm_2^4];

	jacobian=[jac_cen jac_ht jac_wid];
end

splitgauss = splitgauss * ht;
