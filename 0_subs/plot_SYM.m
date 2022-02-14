function plot_SYM(SYM);


debug = 0;

if debug
    disp('startvars :')
    SYM.startvars
    SYM.sim_vars
end

AUX = SYM.AUX;
mask = AUX.mask;
DAT = SYM.DAT;
Nq = SYM.AUX.Nq

for ind = 1:Nq
    if debug; disp(['Q-ind : ', num2str(ind)]); end;
    if sum(mask(:,ind))>0
        qpoint = DAT.HKL_vals(ind,:);
        xdata = DAT.x_dat(mask(:,ind),ind);
        ydata = DAT.y_dat(mask(:,ind),ind);
        edata = DAT.e_dat(mask(:,ind),ind);


        xfit = DAT.x_dat(mask(:,ind),ind);
        yfit = calc_singleQ(AUX, ind);
        yfit = yfit(mask(:,ind))';
        if debug; [ydata, yfit]; end;

        cens = AUX.auxvars(1:end-1, 1, 1);
        hts = AUX.auxvars(1:end-1, 1+ind, 1);   % heights

        hold off;
        errorbar(xdata,ydata,edata,'b--');
        hold on; 
        plot(xfit,yfit,'r-','linewidth',1);     % best-fit line
        plot(cens, hts, 'g*','linewidth',1);    % fitted centers
        plot([0 80],[0 0],'k--');               % x-axis
        
        axis([AUX.eng(1) AUX.eng(end)]);
        title(['column: ',num2str(ind),', q point: ',num2str(qpoint)]);
        xlabel('Energy (meV)');
        ylabel('Intensity (arb. units)');
        legend('Actual Data', 'Fit Data')
        pause
    end
end   % end loop





end % endfunction