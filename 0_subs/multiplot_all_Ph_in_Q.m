function multiplot_all_phonon_single_Q(SYM, ind_q, maxheight);

AUX=SYM.AUX;
SYM=SYM.DAT;
y_fit=calc_DAT_singleQ(AUX,SYM,ind_q);

if ~exist('maxheight'); maxheight=[]; end
centers=AUX.auxvars(:,1,1);
widths=AUX.auxvars(:,1,2)+AUX.auxvars(:,ind_q+1,2);
heights=AUX.auxvars(:,[2:end],1).^2;
asymm=AUX.peak_asymmetry;

xdat=SYM.xdat;
ydat=SYM.ydat;
edat=SYM.edat;

i1=find(AUX.mask,1,'first');
i2=find(AUX.mask,1,'last');
%[x,i1]=min(xdat(AUX.mask));
%[x,i2]=max(xdat(AUX.mask));
%user_range=[xdat(i1,1) xdat(i2,1)];
user_range=[0 6];
eng=xdat(:,1);
eng=eng(:)';		% force eng to row vector


% === plot data and errorbars ===
subplot(4,1,[1 2 3]);
hold off
errorbar(eng,ydat(:,ind_q),edat(:,ind_q),'o')
%titlestr=AUX.good_Q{ind_q};
titlestr=[];
titlestr(regexp(titlestr,'_'))=' ';
title([ num2str(ind_q) ': ' titlestr])
hold on;


% === background (plots straight line between E=0 and point past user range) ===
%user_range=[min(xdat(AUX.mask)) max(xdat(AUX.mask))];
%bg_1=const(ind_q);
%bg_2=const(ind_q) + slope(ind_q)*2*user_range(2);
%plot([0 2*user_range(2)], [bg_1 bg_2], 'k-','linewidth',2);

% scale axes
u_range=0.1*( user_range(2) - user_range(1) );
x_range=[user_range(1)-u_range user_range(2)+u_range];

norm_dat=ydat./edat;
norm_dat=ydat.*(norm_dat>5);
[norm_sort,ind_sort]=sort(conv_nan_to_zero(norm_dat(:)),1,'descend');
%y_range=[-1 1.25*max(max(ydat))];
if isempty(maxheight)
	maxheight=1.25*mean(ydat(ind_sort(1:5)));
end

%disp('flush');fflush(stdout);
y_range=[-1 maxheight];


% === all phonons 
height_vec=heights(:,ind_q);
energy=linspace(eng(1),eng(end),20*length(eng)+1);
all_phonons=zeros(AUX.Nph,length(energy));
w1= widths * asymm/(asymm+1);
w2= widths * 1/(asymm+1);
for ind=1:AUX.Nph
	peak=calc_splitgauss_JAC_fast(energy,centers(ind),1,w1(ind),w2(ind));
	peak=height_vec(ind)*peak(:)';
	all_phonons(ind,:)=peak;
end
plot(energy,all_phonons,'r-','linewidth',1);

% === plot fit ===
%plot(eng, spline(eng,y_fit(:,ind_q), energy), 'r-','linewidth',2);
plot(eng,y_fit,'r-','linewidth',2)
axis([x_range y_range])

% plot diff
subplot(4,1,4); 
hold off;
diffnorm= (ydat(:,ind_q) - y_fit(:,ind_q))./edat(:,ind_q);
plot(eng,diffnorm,'g-','linewidth',2); 
hold on
plot( [eng(1) eng(end)],[0 0],'k-', [eng(1) eng(end)],[2 2],'k--', [eng(1) eng(end)],[-2 -2],'k--', 'linewidth',1);
y_range_diff=[-3 3];
axis([x_range y_range_diff])
set(gca,'xticklabel',[])
%hold off;

