function out=Zezin2015_Na2SO4EOS(input)
% function to evaluate volumetric properties of Na2SO4
% Volumetric Properties of Na2SO4?H2O and Na2SO4?NaCl?H2OSolutions to 523.15 K, 70 MPa
%Denis Zezin, Thomas Driesner, and Carmen Sanchez-Valle
%DOI: 10.1021/je501152a
%J. Chem. Eng. Data 2015, 60, 1181?1192
% Usage: out=Zezin2015_Na2SO4EOS(input) where
%        input is the standard gridded {P,T,m} or scattered [P(:) T(:) m(:)] in MPa K and mol/kg
%        output is structure with appropriate results in appropriate arrays
% calls Fernandez1997dielectric and IAPWS95 to get properties of water
% JMB 2017

% organize input.  Note that P,T,m input arrays are converted into vectors and
% output is reshaped into arrays at the end - This simplifies calculations
% - don't need to track differences between gridded and scattered calculations
    if(iscell(input))
        P=input{1};
        T=input{2};
        m=eps+input{3};
        [pm,tm,mm]=ndgrid(P,T,m);
        nP=length(P);
        nT=length(T);
        nm=length(m);
        tm=tm(:);
        pm=pm(:);
        mm=mm(:);
        npts=length(mm);
        results=Fernandez1997dielectric({P,T});  % get properties of water in cell structure
        Av=results.Av(:)*ones(1,nm);
        Av=Av(:);
        rho=results.rho(:)*ones(1,nm);
        rho=rho(:);
        Vw=1e6*rho.^-1; %factor to convert m^3 to cm^3
    else
        pm=input(:,1);
        tm=input(:,2);
        mm=eps+input(:,3);
        nP=length(pm);
        nT=1;
        nm=1;
        npts=nP;
        results=Fernandez1997dielectric([pm tm]); % get properties of water for scattered data
        Av=results.Av(:);
        rho=results.rho(:);
        Vw=1e6*rho.^-1;  % factor to convert m^3 to cm^3
    end
    
% the physical constants       
R=8.314472;
a1=2;
b=1.2;
nr=1.5;
Po=1;
To=1;
Tr=298.15;
R=8.3144621;
Ms=0.14204;
%Mw=18.0152;
%omega=1000/Mw;
fac=(1+Ms*mm); % factor converting between per kg of water and per kg of solution

I=3*mm; %  I = .5*sum(ni*zi^2)
sqrtI=sqrt(I);
sqrtIr=sqrt(3*1.5);
x=a1*sqrtI;
xr=a1*sqrtIr;
g1 =2*x.^-2.* (1-(1+ x).*exp(-x));
g1r=2*xr.^-2.*(1-(1+xr).*exp(-xr));

tm=tm'; % temporarily transpose tm and pm to make the basis function creation easy to visualize
pm=pm';
basis1=[  % basis functions for Bmx(P,T) and Cmx(P,T)
    ones(1,npts)
    log(tm/Tr)
    (tm-Tr)/To
    To*(620-tm).^-1
    To*(tm-227).^-1
    2*(pm/Po)
    2*(pm/Po).*log(tm/Tr)
    2*(pm/Po).*(tm-Tr)/To
    To*2*(pm/Po).*(620-tm).^-1
    To*2*(pm/Po).*(tm-227).^-1
    ]';
basis2=[  % basis functions f(P,T) for apparent volume at reference concentration 
    1e2*ones(1,npts)
    tm/To
    1e-2*(tm/To).^2
    1e-5*(tm/To).^3
    (pm/Po)
    1e-2*(pm/Po).*(tm/To)
    1e-4*(pm/Po).*(tm/To).^2
    1e-2*(pm/Po).^2
    1e-4*(pm/Po).^2.*(tm/To)
    zeros(1,npts)
    ]';
tm=tm'; % reverse temporary transpose
%pm=pm';

% here are the model parameters from Table 3
parm=[
 6.05954e-4     1.81948e-2    -9.02041e-5     7.85265
-1.33035e-2    -1.22328e-1     2.23418e-3    -9.20087e-1
 1.82630e-5     4.13138e-4     7.54851e-8     2.28824e-1
 2.40232e-1    -4.38235       -1.26843e-1    -7.31316e-2
-6.77293e-2    -1.64437e-1     3.42505e-2    -6.61264e-1
-5.75614e-6    -3.30824e-5    -4.63903e-6     3.58390e-1
 3.61020e-4    -2.45357e-4    -3.60229e-5    -6.66318e-2
-8.08348e-7    -1.11329e-6     6.77063e-8    -1.94421e-1
-1.28047e-3     4.58964e-2     9.04198e-4     6.09792e-2
 5.36298e-4    -8.40841e-3     1.59321e-4      0];

% solve for the Pitzer parameters
beta0=basis1*parm(:,1);
beta1=basis1*parm(:,2);
Cmx=basis1*parm(:,3);
Vovn=basis2*parm(:,4);
Bmx_v=beta0+beta1.*g1;
Bmx_vr=beta0+beta1.*g1r;
V_phi=Vovn-Vw/nr + 6/2/b*Av.*log((1+b*sqrtI)./(1+b*sqrtIr))  + 4*R*tm.*(mm.*Bmx_v-nr*Bmx_vr+2*(mm.^2-nr^2).*Cmx);
V1= 3/b*log(1+b*sqrtI).*Av;
V2=4*R*tm.*(mm.*Bmx_v);
V3=2*4*R*tm.*mm.^2.*Cmx;
Vo=Vovn-Vw/nr + 6/2/b*Av.*log((1+b*sqrtIr).^-1)  + 4*R*tm.*(-nr*Bmx_vr+2*(-nr^2).*Cmx);
dat=(V_phi-V1)./tm/R/4;
Vs=(Vw+V_phi.*mm)./fac;
rhos=1e6*Vs.^-1;  % factor to convert from cm^3 to m^3

out.Vphi=reshape(V_phi,nP,nT,nm);
out.rho=reshape(rhos,nP,nT,nm);
out.Av=reshape(Av,nP,nT,nm);
out.iapws=results.iapws;
out.Vo=reshape(Vo,nP,nT,nm);
out.V1=reshape(V1,nP,nT,nm);
out.V2=reshape(V2,nP,nT,nm);
out.V3=reshape(V3,nP,nT,nm);
out.dat=reshape(dat,nP,nT,nm);

out.beta0=reshape(beta0,nP,nT,nm);
out.beta1=reshape(beta1,nP,nT,nm);
out.Bmx_v=reshape(Bmx_v,nP,nT,nm);
out.c0=reshape(Cmx,nP,nT,nm);
