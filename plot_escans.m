function plot_escans(EXTR)

if ischar(EXTR);
	disp('loading')
	EXTR=load(EXTR);
end

xdat=EXTR.x_dat;
ydat=EXTR.y_dat;
edat=EXTR.e_dat;

ydat(find(ydat==0))=NaN;

qs=size(EXTR.HKL_vals,1)



for ind=1:qs
	errorbar(xdat(:,ind),ydat(:,ind),edat(:,ind),'o');
	hold on; plot([-10 70],[0 0],'k--'); hold off;
	axis([-10 70 -0.2 10])
	pause
	ind
end

