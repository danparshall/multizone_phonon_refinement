function plot_SYM(SYM);



disp('startvars :')
SYM.startvars
SYM.sim_vars

mask = SYM.AUX.mask;
DAT = SYM.DAT;
Nq = SYM.AUX.Nq

for ind = 1:Nq
    if sum(mask(:,ind))>0
        qpoint = DAT.HKL_vals(ind,:);
        xdata = DAT.xdat(mask(:,ind),ind);
        ydata = DAT.ydat(mask(:,ind),ind);
        edata = DAT.edat(mask(:,ind),ind);


        xfit = DAT.xdat(mask(:,ind),ind);
        yfit = calc_JAC_singleQ(SYM.AUX, DAT, ind);
        yfit = yfit(mask(:,ind))';
        [ydata, yfit]


        hold off; errorbar(xdata,ydata,edata,'b--');
        hold on; plot(xfit,yfit,'r-','linewidth',1,[0 80],[0 0],'k--');
        
        axis([DAT.eng(1) DAT.eng(end)]);
        title(['column: ',num2str(ind),', q point: ',num2str(qpoint)]);
        xlabel('Energy (meV)');
        ylabel('Intensity (arb. units)');
        legend('Actual Data', 'Fit Data')
        pause
    end
end   % end loop





end % endfunction