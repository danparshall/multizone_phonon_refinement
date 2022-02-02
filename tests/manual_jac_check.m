function manual_jac_check(SYMS);


for i_sym = 1:length(SYMS)
    SYMS{i_sym}.AUX.freevars(2, :, 1) = 1;
    SYMS{i_sym}.AUX.freevars(2, 1, 2) = 1;



end % sym loop

AUX = SYMS{1}.AUX;
iq = 1;
DELTA = 0.00000001;
%AUX.eng = [0: 0.1, 10];

jac_emp = empirical_jacobian(AUX, iq, DELTA);
[modelout,jac_out] = calc_singleQ(AUX,iq);

size(jac_emp)
size(jac_out)

eng = AUX.eng(AUX.mask(:,iq));
size(eng)

ijac = 8;


%        hold off; errorbar(xdata,ydata,edata,'b--');
        hold on; 
        plot([0 80],[0 0],'k--');
        plot(eng, jac_emp(:, ijac), 'g*','linewidth',1);
        plot(eng, jac_out(:, ijac), 'r-','linewidth',1);
        
        axis([AUX.eng(1) AUX.eng(end)]);
%        title(['column: ',num2str(iq),', q point: ',num2str(qpoint)]);
        xlabel('Energy (meV)');
        ylabel('Intensity (arb. units)');
        legend('Actual Data', 'Fit Data')



%validate_jacobian(SYMS);


end  % endfunc

