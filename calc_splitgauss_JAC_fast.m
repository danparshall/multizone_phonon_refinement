function [splitgauss,jacobian,varargout]=calc_splitgauss_JAC_fast(x_array, cen, ht, hwhm_1, hwhm_2);
% [splitgauss,jacobian]=calc_splitgauss_JAC_fast(x_array, cen, ht, hwhm_1, hwhm_2);
% 	shortest (though confusing) method to calculate splitgauss.
% 	calculates jacobian when nargout==2

x_array=x_array(:);

%splitgauss
splitgauss= ht .* ...
		[exp( -log(2)* ((x_array(x_array<=cen)-cen)/hwhm_1).^2);
		 exp( -log(2)* ((x_array(x_array>cen)-cen)/hwhm_2).^2)];

% jacobian
if nargout==2
	jac_cen= 2*log(2) .* splitgauss .* [(x_array(x_array<=cen)-cen)./hwhm_1.^2; (x_array(x_array>cen)-cen)./hwhm_2.^2];
	jac_ht=splitgauss/ht;
	jac_wid= 2*log(2) .* splitgauss .* [(x_array(x_array<=cen)-cen).^2./hwhm_1.^3; (x_array(x_array>cen)-cen).^2./hwhm_2.^3];

	jacobian=[jac_cen jac_ht jac_wid];
end
