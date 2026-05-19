function results = Bradley1979dielectric(input)
% Usage: results = Bradley1979dielectric(input)
%from Journal of Physical Chemistry, Vol. 83, No. 12, 1979.
%            Daniel J. Bradley and Kenneth S. Pitzer
%           pp.1599-1603
% FORTRAN code by J. CHRISTOPHER PEIPER 29 MARCH, 1982  was consulted and modified
% structure/scatttered input allowed
% Units of MPa, K input in standard {P, T} or [P(:) T(:)] form
% Still a mistqke in the analytic calculation of the Cp DH term (Aj in output) so a
% numerical derivative (Ac) is also returned
% JMB 2017 

% constants from paper:
U1=3.4279e2;
U2=-5.0866e-3;
U3=9.4690e-7;
U4=-2.0525;
U5=3.1159e3;
U6=-1.8289e2;
U7=-8.0325e3;
U8=4.2142e6;
U9=2.1417;

% standard physical constants
%Permittivity of free space,  C^2/J/m
eo=(4e-7*pi*299792458^2)^-1;
%Elementary charge, 
e =1.60217733e-19;
%Boltzmann?s constant, J/K
k= 1.380658e-23;
%Avogadro?s number, 
Na= 6.0221367e23;
R=k*Na;

dT=1e-7;
 
% organize input
    if(iscell(input))
        P=input{1};
        T=input{2};
        [pm,tm]=ndgrid(P,T);
        nP=length(P);
        nT=length(T);
    else
        pm=input(:,1);
        tm=input(:,2);
        nP=length(pm);
        nT=1;
    end
    
 pm=10*pm(:); % factor of 10 to convert P from MPa to bar as used on original paper
 tm=tm(:);
 tm2=tm+dT;

% section calculating the dielectric constant and its derivatives
B=U7 + U8*tm.^-1 +U9*tm ;
BT=U9-U8*tm.^-2 ;
BTT=2.*U8*tm.^-3;

C=U4+U5*(U6+tm).^-1;
CT=-U5*(U6+tm).^-2;
CTT=2.*U5*(U6+tm).^-3 ;

D1K=U1*exp(U2*tm+U3*tm.^2);
D1KT=(U2 + 2.*U3*tm).*D1K ;
D1KTT=(2.*U3+(U2+2.*U3*tm).^2).*D1K ;

D=D1K+C.*log((B+pm)./(B+1000.)) ;
DT=D1KT+C.*((B+pm).^-1 - (B+1000.).^-1).*BT + CT.*log((B+pm).*(B+1000.).^-1) ;
DTT=D1KTT-C.*((B+pm).^-2 - (B+1000.).^-2).*BT.^2+((B+pm).^-1 - (B+1000.).^-1).*(2*BT.*CT + C.*BTT)+CTT.*log((B+pm).*(B+1000.).^-1) ;
DP=10*C.*(B+pm).^-1;
DPP=-1e2*C.*(B+pm).^-2 ;
DTP=-10*C.*BT.*(B+pm).^-2 + CT.*(B+pm).^-1 ;

results.D=reshape(D,nP,nT);
results.DT=reshape(DT,nP,nT);
results.DTT=reshape(DTT,nP,nT);
results.DP=reshape(DP,nP,nT);
results.DPP=reshape(DPP,nP,nT);
results.DTP=reshape(DTP,nP,nT);

% determine second set of dielectric constants for T+dT
B=U7 + U8*tm2.^-1 +U9*tm2 ;
BT=U9-U8*tm2.^-2 ;
C=U4+U5*(U6+tm2).^-1;
CT=-U5*(U6+tm2).^-2;

D1K2=U1*exp(U2*tm2+U3*tm2.^2);
D1KT2=(U2 + 2.*U3*tm2).*D1K2 ;

D2=D1K2+C.*log((B+pm)./(B+1000.)) ;
DT2=D1KT2+C.*((B+pm).^-1 - (B+1000.).^-1).*BT + CT.*log((B+pm).*(B+1000.).^-1) ;

% section calculatiing the limiting slope Debye Huckle terms
% the HAAR EOS for water:
haar= H2O_Haar([pm(:)/10 tm(:)]);
haar2= H2O_Haar([pm(:)/10 tm(:)+dT]);
rho=haar.rho(:);
rho2=haar2.rho(:);
Kt=haar.Kt(:);
alpha=haar.alpha(:);
alpha2=haar2.alpha(:);
dadt=(alpha2-alpha)/dT;

% the standard D-H definitions:   

Aphi=1/3*(2*pi*Na*rho).^.5.*(e.^2*(4*pi*eo*k*D.*tm).^-1).^(3/2);
Agam=Aphi*3;
Av=2*R*tm.*Aphi.*(3e3*DP./D - Kt.^-1)/1e3;
Ah=-6*R*tm.*Aphi.*(1+tm.*DT./D +tm.*alpha/3);

% the equation for Ac here has an error.  I solve numerically below.
% note that the Peiper code leaves out a factor of RT for AH ie the quantity is actually AH/RT
% thus my implimentation 
Ac=Ah.*((2*tm.^-1)-1.5*(DT./D+ tm.^-1 + alpha/3)) - 6*tm.^1.*Aphi.*(DTT./D - (DT./D).^2- tm.^-2 + dadt/3 );

% numerical determination of Ac as derivative of Ah wrt T
% the standard D-H definitions:    
Aphi2=1/3*(2*pi*Na*rho2).^.5.*(e.^2*(4*pi*eo*k*D2.*tm2).^-1).^(3/2);
Ah2=-6*R*tm2.*Aphi2.*(1+tm2.*DT2./D2 +tm2.*alpha2/3);
Ac2=(Ah2-Ah)/dT;

% for reference some of Peiper's unaltered FORTRAN code is reproduced here
%    note that AH is really Ah/RT and AJ is Aj/R in Peiper's code
%       AH=-6.*T*AP*(DT/D +1./T + VT/(3.*V))      
%       AJ=AH*T*(2/T-1.5*(DT/D+1/T+VT/(3*V))) - 6*T**2*AP*(DTT/D -(DT/D)**2-1/T**2+VTT/(3*V)-(VT/V)**2/3)       
%       AV=2.*R*T*AP*(3.*DP/D + VP/V)       
%       AA=AV*(1.+AH/(4.*AP))/T  +2.*R*T*AP*(3.*(DTP/D -DT*DP/D**2)   +VTP/V  -VT*VP/V**2)       
%       AB=-AV**2/(4.*R*T*AP)  +2.*R*T*AP*(3.*(DPP/D -(DP/D)**2) +VPP/V -(VP/V)**2) 

% put results into output structure

results.Aphi=reshape(Aphi,nP,nT);
results.Agam=reshape(Agam,nP,nT);
results.Av=reshape(Av,nP,nT);
results.Ah=reshape(Ah,nP,nT);
results.Aj=reshape(Ac,nP,nT);
results.Ac=reshape(Ac2,nP,nT);
results.haar=haar;
