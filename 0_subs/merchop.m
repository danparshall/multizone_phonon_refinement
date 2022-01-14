function reswid=merchop(ei,frequency,e_trans);
% Dan's brutal modification of MCHOP
% res=chop(ei,frequency);
% includes empirical scale factor of 1.38
scale=1.38;


% Instument details
global x0 xa x1 x2 wa_mm ha_mm wa ha pslit dslat radius rho tjit mod_type s thetam mod_type imod sx sy sz isam gam ia ix idet dd tbin titledata
% first decide if the chopper notation is type and frequency or shorthand
insterr=0;
warning off MATLAB:divideByZero

if length(frequency) > 1 && length(ei) > 1
     disp('easy tiger! set energy range to 0')
     insterr=1;
     return
 end

     %For MERLIN, chop_type=s
chop_type='s';
     x0 = 10.0d0;
     xa = 7.19d0;
     x1 = 1.82d0;
     x2 = 2.5d0; 
     wa_mm = 66.667d0;
     ha_mm = 66.667d0;
     wa = ha_mm / 1000.0d0;
     ha = ha_mm / 1000.0d0;
%     disp('HET setup chosen')
     % chopper details
     if chop_type == 'c'
         chop_par=[ 1.71d0, 0.55d0, 49.0d0,  580.0d0, 0.0d0, 0.0d0 ];
         disp('HET C (100meV) chopper chosen')
         titledata='HET C (100meV)';
     elseif chop_type == 'd'
         chop_par=[ 1.52d0, 0.55d0, 49.0d0,  410.0d0, 0.0d0, 0.0d0 ];
         %chop_par=[.38, 0.02, 10.0d0,  800.0d0, 0.0d0, 0.0d0];
         disp('HET D (50meV) chopper chosen')
         titledata='HET D (50meV)';
     elseif chop_type == 's'
         %chop_par=[2.28d0, 0.55d0, 49.0d0, 1300.0d0, 0.0d0, 0.0d0 ];
         chop_par=[.38, 0.02, 10.0d0,  800.0d0, 0.0d0, 0.0d0];
%         disp('HET S (sloppy) chopper chosen')
%         titledata='HET S (sloppy)';
     elseif chop_type == 'b'
         chop_par= [1.29d0, 0.55d0, 49.0d0,  920.0d0, 0.0d0, 0.0d0 ];
         disp('HET B (200meV) chopper chosen')
         titledata='HET B (200meV)';
     elseif chop_type == 'a'
         chop_par= [0.76d0, 0.55d0, 49.0d0, 1300.0d0, 0.0d0, 0.0d0 ];
         disp('HET A (500meV) chopper chosen')
         titledata='HET A (500meV)';
     else
         disp('Chopper type not recognised')   
     end
     
     % now some moderator details
     % for 300K H2O
     s(1) = 38.6d0;
     s(2) = 0.5226d0;
     s(3) = 0.0d0;
     s(4) = 0.0d0;
     s(5) = 0.0d0;
     th_deg = 26.7d0;
     imod = 2;
     mod_type = 'AP';
     
     % sample details
     sx_mm = 20.0d0;
     sy_mm = 20.0d0;
     sz_mm = 20.0d0;
     isam = 0;
     gam_deg = 0.0d0;
     ia = 0;
     ix = 0;
     
     % detector details
     idet    = 1;
     dd_mm   = 25d0;
     tbin_us = 0.0d0;
     % end of HET parameters 
   


% Convert instrument parameters for the program (set as globals)
pslit  = chop_par(1) ./ 1000.0d0;
dslat  = (chop_par(1) + chop_par(2)) ./ 1000.0d0;
radius = chop_par(3) ./ 1000.0d0;
rho    = chop_par(4) ./ 1000.0d0;
omega = (frequency).*(2.*pi);
tjit   = chop_par(6) * 1.0d-6;


thetam = th_deg*(pi/180.0d0);
% function sigset set a common variable in this case to zero
% sigset (0.0d0, 0.0d0, 0.0d0, 0.0d0)
sx = sx_mm / 1000.0d0;
sy = sy_mm / 1000.0d0;
sz = sz_mm / 1000.0d0;
gam = gam_deg.*pi./180.0d0;

dd = dd_mm ./ 1000.0d0;
tbin = tbin_us .* 1.0d-6;


reswid=get_chop_DAN(ei,omega,e_trans);
reswid=reswid*scale;

