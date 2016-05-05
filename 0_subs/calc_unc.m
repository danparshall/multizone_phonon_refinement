function unc=calc_unc(SYM)
% use


VARS=SYM{1}.VARS;
varsin=VARS.varsin;
ydat=VARS.ydatin;
wdat=VARS.wdatin;		% weights

[func,jac]=calc_model_multiQ(SYM,varsin);

% find Jacobian values that actually change (ignore unused variables)
Jcol=find(sum(jac,1));
jac=jac(:,Jcol);

residuals=func-ydat;

[N_pts,N_vars]=size(jac);
resnorm=residuals.^2 .* wdat;

[Q,R] = qr(jac,0);


%f=inv(full(R));
%covar=sum(f.^2,2);
covar=sum(inv(R).^2,2);

unc=sqrt(covar*sum(resnorm)/(N_pts-N_vars));

