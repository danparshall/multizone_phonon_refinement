function SYMS=make_VARS(SYMS);
% SYMS=make_VARS(SYMS);
% Assembles individual AUX.auxvars into VARS.allvars
%
% auxvars are :
%	(:,:,1)=[cen, height;0, constant BG]
%	(:,:,2)=[wid, res; 0, constant BG];
%
% allvars are: 
%	(:,:,1)=[cen, SYMS1.ht, SYMS2.ht, ...]
%	(:,:,2)=[wid, SYMS1.res, SYMS2.res, ...];
%
% cen/wid are same in both.
%
% col [2:end] of SYMS1.auxvars is 1+[1:AUX.Nq] in allvars
% col [2:end] of SYMS2.auxvars starts after that

ydatin=[];
wdatin=[];
Nvars=[];
genvars=[];
genfree=[];
genLObd=[];
genHIbd=[];
Nqs=[];

for ind=1:length(SYMS)
	AUX=SYMS{ind}.AUX;
	Nqs=[Nqs AUX.Nq];
	DAT=SYMS{ind}.DAT;
	ydatin=[ydatin; DAT.ydat(AUX.mask)];
	wdatin=[wdatin; 1./(DAT.edat(AUX.mask))];
	Nvars=[Nvars length(AUX.indfree)];
	genvars=[genvars, AUX.auxvars(:,[2:end],:)];
	genfree=[genfree, AUX.freevars(:,[2:end],:)];
	genLObd=[genLObd, AUX.bounds_L(:,[2:end],:)];
	genHIbd=[genHIbd, AUX.bounds_H(:,[2:end],:)];
end

VARS.ydatin=ydatin;
VARS.wdatin=wdatin;

VARS.allvars(:,:,1)=[SYMS{1}.AUX.auxvars(:,1,1) genvars(:,:,1)];
VARS.allvars(:,:,2)=[SYMS{1}.AUX.auxvars(:,1,2) genvars(:,:,2)];

VARS.freevars(:,:,1)=[SYMS{1}.AUX.freevars(:,1,1) genfree(:,:,1)];
VARS.freevars(:,:,2)=[SYMS{1}.AUX.freevars(:,1,2) genfree(:,:,2)];

VARS.bndsLO(:,:,1)=[SYMS{1}.AUX.bounds_L(:,1,1) genLObd(:,:,1)];
VARS.bndsLO(:,:,2)=[SYMS{1}.AUX.bounds_L(:,1,2) genLObd(:,:,2)];
VARS.bndsHI(:,:,1)=[SYMS{1}.AUX.bounds_H(:,1,1) genHIbd(:,:,1)];
VARS.bndsHI(:,:,2)=[SYMS{1}.AUX.bounds_H(:,1,2) genHIbd(:,:,2)];


VARS.indfree=find(VARS.freevars);
VARS.Nqs=Nqs;
VARS.Nph=AUX.Nph;

VARS.varsin=VARS.allvars(VARS.indfree);

% connect VARS to first SYMS
SYMS{1}.VARS=VARS;
