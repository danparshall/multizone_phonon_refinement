function plot_SYM(SYM);


debug = 1;

if debug
    disp('startvars :')
%    SYM.startvars

    if isfield(SYM, 'sim_vars')
        disp('sim_vars :')
        SYM.sim_vars
    end
end

AUX = SYM.AUX;
mask = AUX.mask;
DAT = SYM.DAT;
Nq = SYM.AUX.Nq

for ind = 1:Nq
    if debug; disp(['Q-ind : ', num2str(ind)]); end;
    if sum(mask(:,ind))>0
        qpoint = DAT.HKL_vals(ind,:);
%        xdata = DAT.x_dat(mask(:,ind),ind);
        xdata = AUX.eng(mask(:,ind));
        ydata = DAT.y_dat(mask(:,ind),ind);
        edata = DAT.e_dat(mask(:,ind),ind);

%eng_inds = find(mask(:,ind))
%size(eng_inds)
        xfit = AUX.eng(mask(:,ind));
%        xfit = DAT.x_dat(:,ind);
        yfit = calc_singleQ(AUX, ind)';
%        yfit = nan*size(AUX.eng);
%        yfit(mask(:,ind)) = ycalc;
        if debug; [ydata, yfit]; end;

        cens = AUX.auxvars(1:end-1, 1, 1);
        hts = AUX.auxvars(1:end-1, 1+ind, 1);   % heights

        hold off;
        errorbar(xdata,ydata,edata,'b--');
        hold on; 
        plot(xfit,yfit,'r-','linewidth',1);     % best-fit line
        plot(cens, hts, 'g*','linewidth',1);    % fitted centers
        plot([0 80],[0 0],'k--');               % x-axis
        
        axis([AUX.eng(1) AUX.eng(end) -1 10]);
        title(['column: ',num2str(ind),',  q point: [', num2str(qpoint(1)) ', ' num2str(qpoint(2)) ', ' num2str(qpoint(3)) ']']);
        xlabel('Energy (meV)');
        ylabel('Intensity (arb. units)');
        legend('Measured Data', 'Fit Result')
        pause
    end
end   % end loop





end % endfunction