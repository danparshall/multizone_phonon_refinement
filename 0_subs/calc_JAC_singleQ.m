function [modelout,jacout,varargout]=calc_JAC_singleQ(AUX,DAT,idx);
% [modelout,jacout,varargout]=calc_JAC_singleQ(AUX,DAT,idx);
% 	Calculates the data for a single Q.  
%	If requested, calculates and delivers Jacobian as well.
%	option in here to constrain heights.
%	Cycles through, calculates profile and jacobian for each phonon.

if ~exist('idx'); idx=1; end

eng = AUX.eng;
eng=eng(:)';

cen=AUX.auxvars(1:end-1,1,1);
ht=AUX.auxvars(1:end-1,idx+1,1);

%disp([' Zero-peaks: ', num2str(length(find(ht==0)))]);

%%%% WIDTHS
FWHM = AUX.auxvars(1:end-1,1,2) + AUX.auxvars(1:end-1,idx+1,2); % add phonon+resolution width
asymm=AUX.peak_asymmetry;
w1= FWHM * asymm/(asymm+1);
w2= FWHM * 1/(asymm+1);

% Background
constant = AUX.auxvars(end,idx+1,1);
slope = AUX.auxvars(end,idx+1,2);

%%% INIT
modelout=zeros(size(eng));

if nargout==2;		% last two cols are derivatives of background
	jacout=[zeros(length(eng),3*AUX.Nph) ones(length(eng),1) eng(:)];
end


for ind=1:AUX.Nph
	if nargout==2	% when asked for jacobian
		[splitgauss,splitjac]=calc_splitgauss_JAC_fast(eng,cen(ind),ht(ind),w1(ind),w2(ind));

		%disp(['jacsize = ', num2str(size(splitjac))]);
		modelout=modelout+splitgauss';
		splitjac=splitjac(:,[1 2 3]);	% cen ht wid (both widths combined)
		ind_j= (3*ind)-2;
		jacout(:,[ind_j:ind_j+2])=splitjac;

	else			% no jacobian
		splitgauss=calc_splitgauss_JAC_fast(eng,cen(ind),ht(ind),w1(ind),w2(ind));
		modelout=modelout+splitgauss';
	end
end


% tack on the backgrounds to the convolved splitgausses
modelout = modelout + constant;% constant BG
modelout = modelout + (slope*eng);% linear BG



if 1
	disp(['modelout = ', num2str(size(modelout))]);
	disp(['jacout = ', num2str(size(jacout))]);
	validEng = AUX.mask(:,idx);
	disp(['eng = ', num2str(size(eng))]);
	disp(['ydat= ', num2str(size(DAT.ydat(AUX.mask(:,idx),idx)))]);
	disp(['edat= ', num2str(size(DAT.edat(AUX.mask(:,idx),idx)))]);
%DAT.edat(AUX.mask(:,idx),idx)
disp('FLUSH')
fflush(stdout);
pause
	hold off; errorbar(AUX.eng(validEng),DAT.ydat(validEng,idx),DAT.edat(validEng,idx),'b--')
	hold on; plot(eng,modelout,'r-','linewidth',1,[eng(1) eng(end)],[0 0],'k--')
	vec = axis
	pause
end
%disp(['modelout = ', num2str(size(modelout))]);
%disp(['jacout = ', num2str(size(jacout))]);
if system_octave; fflush(stdout); end

