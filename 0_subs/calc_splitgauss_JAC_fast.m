function [splitgauss,jacobian,varargout]=calc_splitgauss_JAC_fast(eng, cen, ht, hwhm_1, hwhm_2);
% [splitgauss,jacobian]=calc_splitgauss_JAC_fast(eng, cen, ht, hwhm_1, hwhm_2);
% 	shortest (though confusing) method to calculate splitgauss.
% 	calculates jacobian when nargout==2
%	splitgauss size is [ len(eng) ]
%	jacobian size is [ len(eng) x 3 ] (because the two widths are combined)

eng=eng(:);

%splitgauss
splitgauss= [exp( -log(2)* ((eng(eng<=cen)-cen)/hwhm_1).^2);
		 	 exp( -log(2)* ((eng(eng>cen)-cen)/hwhm_2).^2)];

% jacobian
if nargout==2
	jac_cen= ht*2*log(2) * splitgauss .* [(eng(eng<=cen)-cen)./hwhm_1.^2; (eng(eng>cen)-cen)./hwhm_2.^2];
	jac_ht = splitgauss;
	jac_wid= ht*2*log(2) * splitgauss .* [(eng(eng<=cen)-cen).^2./hwhm_1.^3; (eng(eng>cen)-cen).^2./hwhm_2.^3];

	jacobian=[jac_cen jac_ht jac_wid];
end

splitgauss = splitgauss * ht;
