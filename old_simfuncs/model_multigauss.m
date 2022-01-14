function [modelout,jacout]=model_multigauss(gdata,eng)
% each row of gdata is center, height, W1, W2

eng=eng(:)';

modelout=zeros(size(eng));
jacout=zeros(length(eng),numel(gdata));

for ind=1:size(gdata,1);
	[splitgauss,splitjac]=calc_splitgauss_JAC_fast(eng,gdata(ind,1),gdata(ind,2),gdata(ind,3), gdata(ind,4));

	modelout=modelout+splitgauss';
	ind_j= (4*ind)-3;
	jacout(:,[ind_j:ind_j+2])=splitjac;
end
%jacout
