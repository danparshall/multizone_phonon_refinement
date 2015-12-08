function [modout,jacout,varargout]=calc_JAC_singleQ(AUX,SYM,idx);
% [modout,jacout,varargout]=calc_JAC_singleQ(AUX,SYM,idx);
% 	Calculates the data for a single Q.  
%	If requested, calculates and delivers Jacobian as well.
%	option in here to constrain heights.

if ~exist('idx'); idx=1; end

eng=SYM.xdat(AUX.mask(:,idx),idx);	%	pulls just valid points
eng=eng(:)';

cen=AUX.auxvars(1:end-1,1,1);
ht=AUX.auxvars(1:end-1,idx+1,1);

%%%% WIDTHS
FWHM = AUX.auxvars(1:end-1,1,2) + AUX.auxvars(1:end-1,idx+1,2); % add phonon+resolution width
asymm=AUX.peak_asymmetry;
w1= FWHM * asymm/(asymm+1);
w2= FWHM * 1/(asymm+1);

% Background
constant = AUX.auxvars(end,idx+1,1);
slope = AUX.auxvars(end,idx+1,2);

%%% INIT
modout=zeros(size(eng));
if nargout==2; jacout=zeros(length(eng),3*AUX.Nph); end

if 1 % normal mode, ht not constrained
	for ind=1:AUX.Nph
		if nargout==2	% when asked for jacobian
			[splitgauss,splitjac]=calc_splitgauss_JAC_fast(eng,cen(ind),ht(ind),w1(ind),w2(ind));
			modout=modout+splitgauss';
			splitjac=splitjac(:,[1 2 3]);					% cen ht wid
			ind_j= (3*ind)-2;
			jacout(:,[ind_j:ind_j+2])=splitjac;
		else			% no jacobian
			splitgauss=calc_splitgauss_JAC_fast(eng,cen(ind),ht(ind),w1(ind),w2(ind));
			modout=modout+splitgauss';
		end
	end

else	% heights constrained to be positive
	for ind=1:AUX.Nph
		if nargout==2	% when asked for jacobian
			[splitgauss,splitjac]=calc_splitgauss_JAC_fast(eng,cen(ind),1,w1(ind),w2(ind));
			modout=modout+splitgauss'*ht(ind)^2;
			splitjac=splitjac*ht(ind)^2;					% scale all by height^2
			splitjac(:,2)=splitjac(:,2)*2/ht(ind);			% scale htjac by 2/ht (so overall, 2*ht)
			ind_j= (3*ind)-2;
			jacout(:,[ind_j:ind_j+2])=splitjac;
		else			% no jacobian
			splitgauss=calc_splitgauss_JAC_fast(eng,cen(ind),1,w1(ind),w2(ind));
			modout=modout+splitgauss'*ht(ind)^2;
		end
	end
end

%tack on the backgrounds to the convolved splitgausses
modout = modout + constant;% constant BG
modout = modout + (slope*eng);% linear BG

if nargout==2 % add derivatives of background to the jacobian
	jacout = [jacout ones(length(eng),1) eng'];
end

if 0
	hold off; errorbar(eng,SYM.ydat(AUX.mask(:,idx),idx),SYM.edat(AUX.mask(:,idx),idx),'b--')
	hold on; plot(eng,modout,'r-','linewidth',1)
	pause
end

fflush(stdout);
