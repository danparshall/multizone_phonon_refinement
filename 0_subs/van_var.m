function [v_van,v_van_m,v_van_ch,v_van_jit,v_van_ya,v_van_x,v_van_y,v_van_xy,v_van_dd]=van_var(ei,omega,etrans)
global x0 xa x1 x2 wa_mm ha_mm wa ha pslit dslat radius rho tjit mod_type s thetam mod_type imod sx sy sz isam gam ia ix idet dd tbin
warning off MATLAB:divideByZero
% VAN_VAR (x0, xa, x1, x2, wa, s, thetam, imod, pslit, radius, rho, omega, tjit,
%      &                      sx, sy, sz, isam, ix, ia, dd, idet, ei, eps, phi, gam, v_van)
% !
% !  Calculates the inelastic vanadium peak variance seen in the detector including the effects of rotor jitter.
% !  DOES NOT INCLUDE EFFECT OF TIME_BIN WIDTH.
% !  All distances in metres; WA full width. All times are in S. THETAM in rad, OMEGA in rad/S.
% !  The individual components of the width are returned in the common block /VAN_VARS/ by the
% ! routine VAN_CALC.
% !

%! Parameters:
nsp = 10;
% ! Internal variables:
% 	integer ifail
% 	double precision sp(nsp), wi, wf, tsqmod, tsqchp, tsqjit, atten, dx, dy, dz, v_x, v_y, v_z, v_xy,
%      &                 deld, sigd, sigdz, sigdd, effic, v_dd

% eps energy transfer =0 
%!----------------------------------------------------------------------------------------------------------
if nargin == 3
    eps=etrans;	
    
else
    eps=0;
end
 
       
wi = 0.69468875.*sqrt(ei);
wf = 0.69468875.*sqrt(ei-eps);

%    
% !  get a load of widths:
% !  ---------------------
% !  moderator:
if(imod == 0) 
    tsqmod=tchi(s(1)./1000.0d0, ei);
else if(imod == 1)
        tsqmod=tikeda(s(1), s(2), s(3), s(4), s(5), ei);
    else if(imod == 2)
            tsqmod=tchi_2 ((s(1)./1000.0d0), (s(2)./1000.0d0), ei);
            
        end
    end
end

    
%!  chopper:
[tsqchp,ifail]=tchop(omega, ei);
ifail;
if (ifail ~= 0) ;
    tsqchp = 0.0d0;
end

%!  chopper jitter:
tsqjit = tjit.^2;

%!  sample: it appears ia is always set to zero!
% if(ia =~ 0)
% sam4 (sx, sy, sz, isam, ix, sp, nsp, wi, wf, gam, phi, 6, ifail,atten, dx, dy, dz, v_x, v_y, v_z, v_xy)
% else
[v_x, v_y, v_z]=sam0(6);
atten=1.0d0;
dx=0.0d0;
dy=0.0d0;
dz=0.0d0;
v_xy=0.0d0;

% !  detector:
[deld, sigd, sigdz, sigdd, effic]=detect2(1.0d0, 1.0d0, wf, 6);
v_dd=sigdd.^2;

% !  now calculate the inelastic vanadium width:
% !  -------------------------------------------
phi=10;
[v_van,v_van_m,v_van_ch,v_van_jit,v_van_ya,v_van_x,v_van_y,v_van_xy,v_van_dd]=van_calc(tsqmod, tsqchp, tsqjit, v_x, v_y, v_xy, v_dd, ei, eps, phi,omega);

%moderator functions    
    
% CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
% C
% C     SUBROUTINE TIKEDA(S1,S2,B1,B2,EMOD,EI,TAUSQR)
% C     =============================================
% C  Works out the variance of the moderator pulse for a moderator
% C  described by an Ikeda-Carpenter time pulse. Answer in sec**2.
% C
% C     S1  constant in expression for macroscopic xsect        R*8
% C        of moderator. Units: m**-1
% C     S2  gradient of wavelength**2 in expression for xsect   R*8
% C        Units: (m*Ang)**-1
% C     B1  inverse time constant for moderator storage term    R*8
% C        in the range E<130 meV. Units: mms**-1
% C     B2  inverse time constant for E>130 meV                 R*8
% C        Units: mms**-1
% C   EMOD  swap over energy from chi squared to storage        R*8
% C        term (meV)
% C     EI  energy of neutrons (meV)                            R*8
% C TAUSQR  value of the variance of Ikeda/Carpenter            R*8
% C        function at the specified energy. Returned
% C        in sec**2
% C
% C ENTRY
% C =====
% C   All must be provided, except for TAUSQR which is supplied
% C  by the the routine.
% C
% C EXIT
% C ====
% C   All unchanged except TAUSQR which is now assigned.
% C
% CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
function tausqr=tikeda(S1,S2,B1,B2,EMOD,ei)
global x0 xa x1 x2 wa_mm ha_mm wa ha pslit dslat radius rho tjit mod_type s thetam mod_type imod sx sy sz isam gam ia ix idet dd tbin

SIG=sqrt( (S1.*S1) + ((S2.*S2.*81.8048)./ei) );
A = 4.37392D-4 .* SIG .* sqrt(ei);
for j=1:length(ei)
    if (ei(j) > 130.0);
        B(j)=B2;
    else
        B(j)=B1;
    end
   

R=exp(-ei./EMOD);
tausqr(j)=(3.0./(A.*A)) + (R.*(2.0-R))./(B(j).*B(j));
end
% variance currently in mms**2. Convert to sec**2
tausqr=tausqr.*1.0D-12;



% CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
% C
% C     SUBROUTINE TCHI(DELTA,EI,TAUSQR)
% C
% C  Works out the variance of the moderator pulse for a
% C  Chi**2 function. Answer in S**2.
% C
% C     DELTA    characteristic thickness of moderator. It        R*8
% C             is the distance a neutron of energy EI covers
% C             in a time equal to the FWHH of
% C             the moderator time pulse. It is a constant for
% C             the Chi**2 function. Is 28mm typicaly. DELTA
% C             is needed in metres.
% C     EI       energy of neutrons (meV)                         R*8
% C     TAUSQR   value of the variance of Chi**2 function         R*8
% C             at the specified energy. Returned in sec**2
% C
% C ENTRY
% C =====
% C  DELTA,EI are needed.
% C  TAUSQR set by the program
% C
% C EXIT
% C ====
% C  DELTA,EI unchanged.
% C  TAUSQR assigned.
% C
% CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

function [tausqr]=tchi(DELTA,ei)   
global x0 xa x1 x2 wa_mm ha_mm wa ha pslit dslat radius rho tjit mod_type s thetam mod_type imod sx sy sz isam gam ia ix idet dd tbin

VEL=437.392.*sqrt(ei);
tausqr=( (DELTA./1.96)./ VEL ).^2;


% CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCccccccc
% C
% C     SUBROUTINE TCHI_2(DELTA_0,DELTA_G,EI,TAUSQR)
% C
% C  Works out the variance of the moderator pulse for a
% C  Chi**2 function with an energy dependant width of the form:
% C
% C     DELTA  =  DELTA_0  +  DELTA_G * SQRT(EmeV)
% C
% C   where:
% C
% C     DELTA    characteristic thickness of moderator. It
% C             is the distance a neutron of energy EI covers
% C             in a time equal to the FWHH of
% C             the moderator time pulse. It is a constant for
% C             the Chi**2 function. Is 28mm typicaly.
% C
% C  Answer in S**2.
% C
% C     DELTA_0  in metres                                        R*8
% C     DELTA_G  in metres/(meV)**0.5                             R*8
% C     EI       energy of neutrons (meV)                         R*8
% C     TAUSQR   value of the variance of Chi**2 function         R*8
% C             at the specified energy. Returned in sec**2
% C
% C ENTRY
% C =====
% C  DELTA_0,DELTA_G,EI are needed.
% C  TAUSQR set by the program
% C
% C EXIT
% C ====
% C  DELTA_0,DELTA_G,EI unchanged.
% C  TAUSQR assigned.
% C
% CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

function [tausqr]=tchi_2(DELTA_0,DELTA_G,ei)
global x0 xa x1 x2 wa_mm ha_mm wa ha pslit dslat radius rho tjit mod_type s thetam mod_type imod sx sy sz isam gam ia ix idet dd tbin

VEL=437.392.*sqrt(ei);
tausqr=( ( (DELTA_0+DELTA_G.*sqrt(ei))./1.96) ./ VEL ).^2;

% end of moderator functions

function [tausqr,ierr]=tchop(omega,ei)
global x0 xa x1 x2 wa_mm ha_mm wa ha pslit dslat radius rho tjit mod_type s thetam mod_type imod sx sy sz isam gam ia ix idet dd tbin

% !
% !   Calculates the variance of the time pulse through a rotor at
% !  any energy. Answer in seconds**2
% !   New version 15/1/90
% !
% !    p       slit thickness (m)                              R*8
% !    R       slit package diameter (m)                       R*8
% !    rho     slit radius of curvature (m)                    R*8
% !    w       angular frequency of rotor (rad/sec)            R*8
% !    ei      energy the rotor has been phased for (meV)      R*8
% !    tausqr  variance of time width of rotor                 R*8
% !    ierr    error indicator                                 integer
% !               =0  no problems
% !               =1  if no transmission;TAUSQR set to zero
% !
% !  Entry
% !  =====
% !    p,R,rho,w,ei must be supplied
% !    tausqr,ierr assigned by the program, so entry values are
% !   unimportant.
% !
% !  Exit
% !  ====
% !    p,R,rho,w,ei unchanged
% !    tausqr,ierr assigned by the program
% !
p=pslit;  
R=radius; 
rho=rho  ; 
w=omega  ;
ei=ei;


if (p == 0.0d0 && R == 0.0d0 & rho == 0.0d0) 
    ierr=1;
    tausqr = 0.0d0;
end

% !  Calculate parameter gam:
% ! --------------------------
      veloc=437.392d0.*sqrt(ei);
      gammm=( 2.0d0.*(R.^2)./p ) .* abs(1.0d0./rho - 2.0d0.*w./veloc);

%!  Find regime and calculate variance:
%! -------------------------------------
for j=1:length(ei)
    groot=0;
    if (gammm(j) >= 4.0d0)
        ierr=1;
        tausqr=0.0d0;
    else
        ierr=0;
        if (gammm(j) <= 1.0d0)
            gsqr(j)=(1.0d0-(gammm(j).^2).^2 ./10.0d0) ./ (1.0d0-(gammm(j).^2)./6.0d0);
        else
            groot=sqrt(gammm(j));
            gsqr(j)=0.6d0.*gammm(j).*((groot-2.0d0).^2).*(groot+8.0d0)/(groot+4.0d0);
        end
        tausqr(j)=( (p./(2.0d0.*R.*w)).^2 ./ 6.0d0) .* gsqr(j);
    end
end



% ! SAM0: Given the sample type the routine calculates the 
% !      variances given the characteristic sample dimensions.
% !      The sample types are all non-absorbing.
% !
% !     SX     x-axis charcteristic dimension             R*8
% !     SY     y-axis charcteristic dimension             R*8
% !     SZ     z-axis charcteristic dimension             R*8
% !     ISAM   sample type : 0 non-absorbing plate        integer
% !                           (SX... full widths )
% !                          1 non-absorbing ellipsoid
% !                           (SX... FULL widths)
% !                          2 hollow cylinder
% !                           (SX = external radius
% !                            SY = internal radius
% !                            SZ = full height)
% !                          3 sphere
% !                           (SX= radius
% !                            SY, SZ ignored)
% !                          4 solid cylinder
% !                           (SX = external radius
% !                            SY = ignored
% !                            SZ = full height)
% !     IOUT   output stream number                       integer
% !     VARX   variance for x-axis                        R*8
% !     VARY   variance for y-axis                        R*8
% !     VARZ   variance for z-axis                        R*8
% !     IERR   =0 no problems                             integer
% !             1 invalid sample type
% !
% ! ENTRY
% ! =====
% !   SX,SY,SZ,ISAM,IOUT must be provided
% !
% ! EXIT
% ! ====
% !   SX,SY,SZ,ISAM,IOUT unchanged
% !   VARX,VARY,VARZ,IERR returned
% !
function [varx, vary, varz]=sam0(iout)
global x0 xa x1 x2 wa_mm ha_mm wa ha pslit dslat radius rho tjit mod_type s thetam mod_type imod sx sy sz isam gam ia ix idet dd tbin


ierr = 0;
if (isam == 0)         % ! plate sample
    varx = sx.^2 ./ 12.0d0;
    vary = sy.^2 ./ 12.0d0;
    varz = sz.^2 ./ 12.0d0;
else if (isam == 1)     % ! ellipsoidal sample
        varx = sx.^2 ./ 20.0d0;
        vary = sy.^2 ./ 20.0d0;
        varz = sz.^2 ./ 20.0d0;
    else if (isam == 2)    % ! ;hollow cylinder
            varx = (sx.^2 + sy.^2) ./ 4.0d0;
            vary = (sx.^2 + sy.^2) ./ 4.0d0;
            varz = sz.^2 ./ 12.0d0;
        else if (isam == 1)    % ! spherical sample
                varx = sx.^2 ./ 20.0d0;
                vary = sx.^2 ./ 20.0d0;
                varz = sx.^2 ./ 20.0d0;
            else if (isam == 4)    %  ! solid cylinder
                    varx = sx.^2 ./ 4.0d0;
                    vary = sx.^2 ./ 4.0d0;
                    varz = sz.^2 ./ 12.0d0;
                end
            end
        end
    end
end
 
% ! sub DETECT2
% !  Entry:
% !     wd  width of detector (m)
% !     hd  height of detector (m)
% !     dd  depth of deitector (m)
% !    idet detector type
% !      =0  He gas detector
% !            dd     is taken as the DIAMETER
% !            wd/dd  gives the number binned sideways (nearest integer)
% !
% !      =1  He gas detector
% !            As above except that wd is treated as the width of a hat
% !            function. This isuseful if the gas tubes are not spaced
% !            as their diameters
% !
% !      =2  scintillator of the Davidson type
% !
% !     wf  final wave-vector (Ang-1)
% !   iout  output stream for error messages
% !
% !  Exit:
% !   delta   shift in nominal position of detector from the midway
% !          position (e.g. gas-tube from centre, scintillator from
% !          dd/2)   (m)
% !    sigd   variance in the width   (m)
% !   sigdz   variance in the height  (m)
% !   sigdd   variance in the depth   (m)
% !   effic   efficiency (0<effic<1)
% !    ierr  error flag
% !      =0  no problems
% !      =1  problems
% !
% !
function [delta,sigd,sigdz,sigdd,effic]=detect2(wd,hd,wf,iout)
global x0 xa x1 x2 wa_mm ha_mm wa ha pslit dslat radius rho tjit mod_type s thetam mod_type imod sx sy sz isam gam ia ix idet dd tbin

unity=1.0;
ierr=0;
% seems this is the only useable function
% !  He cylindrical detectors binned together
% ! -----------------------------------------
rad=dd./2.0;
atms=10.0;
t2rad=0.063;

%don't call this function at the mo and approximate some values for it
%[effic,delta,ddsqr,v_dd,v_d]=detect_he(wf,rad,atms,t2rad)
effic=.5;
delta=0;
ddsqr=0.25;
v_dd=0.01;
v_d=0.01;
    sigd =wd./sqrt(12.0);
    sigdz=hd./sqrt(12.0);
    sigdd=sqrt(v_dd);
	
% !  function van_calc
% !  Calculates the inelastic vanadium width seen in the detector including the effects of rotor jitter.
% !  DOES NOT INCLUDE EFFECT OF TIME_BIN WIDTH.                                                                 !  All distances in metres; WA full width. All times are variances in S. THETAM in rad, OMEGA in rad/S.
% !  The individual components of the width are returned in the common block /VAN_VARS/
% !
function [v_van,v_van_m,v_van_ch,v_van_jit,v_van_ya,v_van_x,v_van_y,v_van_xy,v_van_dd]=van_calc(v_mod, v_ch, v_jit, v_x, v_y, v_xy, v_dd,ei, eps, phi,omega)
global x0 xa x1 x2 wa_mm ha_mm wa ha pslit dslat radius rho tjit mod_type s thetam mod_type imod sx sy sz isam gam ia ix idet dd tbin

veli=437.3916d0.*sqrt( ei );
velf=437.3916d0.*sqrt( ei-eps );
rat=(veli./velf).^3;


tanthm=tan(thetam);

am   = -(x1+rat.*x2)./x0;
ach  = (1.0d0 + (x1+rat.*x2)./x0);
g1 = (1.0d0 - (omega.*(x0+x1).*tanthm./veli)); 
g2 = (1.0d0 - (omega.*(x0-xa).*tanthm./veli) ) ;
f1 =  1.0d0 + ((x1./x0).*g1);
f2 =  1.0d0 + ((x1./x0).*g2);
gg1 = g1 ./ ( omega.*(xa+x1) ) ;
gg2 = g2 ./ ( omega.*(xa+x1) ) ;
ff1 = f1 ./ ( omega.*(xa+x1) ) ;
ff2 = f2 ./ ( omega.*(xa+x1) ) ;
aa = ( (cos(gam)./veli) - (cos(gam-phi)./velf) ) - (ff2.*sin(gam));
bb = ((-sin(gam)./veli) + (sin(gam-phi)./velf) ) - (ff2.*cos(gam));
aya  = ff1 + ((rat.*x2./x0).*gg1);
ax   = aa  - ((rat.*x2./x0).*gg2.*sin(gam));
ay   = bb  - ((rat.*x2./x0).*gg2.*cos(gam));
a_dd  = 1.0d0./velf;

v_van_m  = am.^2  .* v_mod;
v_van_ch = ach.^2 .* v_ch;
v_van_jit= ach.^2 .* v_jit;
v_van_ya = aya.^2 .* (wa.^2./12.0d0);
v_van_x  = ax.^2  .* v_x;
v_van_y  = ay.^2  .* v_y;
v_van_xy = ax.*ay  .* v_xy;
v_van_dd = a_dd.^2.* v_dd;

v_van = (v_van_m + v_van_ch + v_van_jit + v_van_ya);%    + v_van_x);% + v_van_y);%  + v_van_xy);%  + v_van_dd);


