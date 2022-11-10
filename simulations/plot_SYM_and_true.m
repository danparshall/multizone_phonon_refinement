function plot_SYM_and_true(SYM, true_y);


debug = 0;

if debug
    disp('startvars :')
    SYM.startvars
    disp('sim_vars :')
    SYM.sim_vars
end

AUX = SYM.AUX;
mask = AUX.mask;
DAT = SYM.DAT;
Nq = SYM.AUX.Nq;

if debug; disp(['Mask size : ', num2str(size(mask))]); end;

for ind = 1:Nq
    if debug; disp(['Q-ind : ', num2str(ind)]); fflush(stdout); end;
    if sum(mask(:,ind))>0
        qpoint = DAT.Q_hkl(ind,:);
        xdata = DAT.x_dat(mask(:,ind),ind);
        ydata = DAT.y_dat(mask(:,ind),ind);
        edata = DAT.e_dat(mask(:,ind),ind);


        xfit = DAT.x_dat(mask(:,ind),ind);
        yfit = calc_singleQ(AUX, ind);
%        size(yfit)
%        yfit = yfit(mask(:,ind))';
        if debug; [ydata, yfit(:)]; end;

        cens = AUX.auxvars(1:end-1, 1, 1);
        hts = AUX.auxvars(1:end-1, 1+ind, 1);   % heights

        hold off;
        errorbar(xdata,ydata,edata,'b--');
        hold on; 
        plot(xfit,yfit,'r-','linewidth',1);     % best-fit line

        if exist('true_y')
            heights = SYM.sim_vars(:, 3:end);
%            plot(SYM.sim_vars(:,1), SYM.sim_vars(:, 2+ind), 'g*');   % true underlying peak centers/heights
            plot(SYM.sim_vars(:,1), heights(:,ind), 'g*');   % true underlying peak centers/heights
            plot(xfit, true_y(mask(:,ind), ind), 'g:', 'linewidth', 2);  % ideal scan, given the above
        end

        range_xmax = max([0.7*max(AUX.eng), 1.2*max(cens)]);
        plot([0 range_xmax],[0 0],'k--');               % x-axis
        plot(cens, hts, 'r*','linewidth',1);    % fitted centers
        
        axis([DAT.eng(1) range_xmax]);
        title(['column: ',num2str(ind),', q point: ',num2str(qpoint)]);
        xlabel('Energy (meV)');
        ylabel('Intensity (arb. units)');
        legend('Noisy Data', 'Fit Data', 'Ground Truth')
        pause
    end
end   % end loop





end % endfunction