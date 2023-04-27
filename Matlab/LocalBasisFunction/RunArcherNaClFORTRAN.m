function Archer=RunArcherNaClFORTRAN(input)
% this function is a frontend to the batch FORTRAN code of Archer. The
% input conditions of P (MPa), T (K), and m (mol/kg) are in a cell. which is then converted
% into a list of PTm points.  
% Usage:
%     Output=RunArcherNaClFORTRAN(PTm)
%   Output=  struct with fields:
% 
%        gam: dimesionless
%        phi: dimesionless
%         aw: dimesionless
%         AV: cc/mol
%      alpha: K^-1
%       beta: 
%         Hr: J/kg
%        ACp: J/mol/K
%       rhos: kg/m^3
%       rhow: kg/m^3
%     alphas: K^-1
%     alphaw: K^-1
%      betas: 
%      betaw: 
%        Cps: J/kg/K
%        Cpw: J/kg/K
%         Go: J/mol
%         Ho: J/mol
%         So: J/mol/K
%       Aphi: dimesionless
%         Av: dimesionless
%         Ac: dimesionless
%         Ah: dimesionless
%       diel: 77.718
%         Gw: J/kg
%         Sw: J/kg/K
%          G: J/kg
%         Gs: J/kg
% JMB 2017 

[pm,tm,mm]=ndgrid(input{1},input{2},input{3});
Pflg=0;
Tflg=0;
mflg=0;
nP=length(input{1});
nT=length(input{2});
nm=length(input{3});
if nP==1, Pflg=1;end
if nT==1, Tflg=1;end
if nm==1, mflg=1;end

PTm=[pm(:) tm(:) mm(:)];
R=8.3144;
omega=1000/18.0152;
MW=.058443;
Goo=-9.045;

output=Archer_NACL(PTm);

% 1 T, K    
% 2 P, MPa   
% 3 M(mol. kg-1)  
% 4 sol.act.coeff.  
% 5 solv.act.coeff.  
% 6 solv.activity  
% 7 V/(cm3.mol-1)  
% 8 expans./(cm3.mol-1.K-1)  
% 9 compr./(cm3.mol-1.MPa-1)  
% 10 rel. enthalpy/(kJ.mol-1)  
% 11 Cp/(J.K-1.mol-1)  
% 12 Soln dens/(g.cm-3)  
% 13 Solv dens/(g.cm-3)  
% 14 Soln expans/(K-1)  
% 15 Solv expans/(K-1)  
% 16 Soln compres/(MPa-1)  
% 17 Solv compres/(MPa-1)  
% 18 Soln Cp/(J.K-1.mol-1)  
% 19 Solv Cp/(J.K-1.mol-1)  
% 20 G0-G0(Tr,pr)/kJ.mol-1  
% 21 H0-H0(Tr,pr)/kJ.mol-1
% 22 S0-S0(Tr,pr)/J.K-1.mol-1
% 23 Aphi
% 24 D
% 25 Gw
%26 Sw
Archer.gam=squeeze(squeeze(squeeze(reshape(output(:,1),nP,nT,nm))));
Archer.phi=squeeze(squeeze(squeeze(reshape(output(:,2),nP,nT,nm))));
Archer.aw=squeeze(squeeze(squeeze(reshape(output(:,3),nP,nT,nm))));
Archer.AV=squeeze(squeeze(squeeze(reshape(output(:,4),nP,nT,nm))));
Archer.alpha=squeeze(squeeze(squeeze(reshape(output(:,5),nP,nT,nm))));
Archer.beta=squeeze(squeeze(squeeze(reshape(output(:,6),nP,nT,nm))));
Archer.Hr=1e3*squeeze(squeeze(squeeze(reshape(output(:,7),nP,nT,nm))));
Archer.ACp=1e3*squeeze(squeeze(squeeze(reshape(output(:,8),nP,nT,nm))));
Archer.rhos=1e3*squeeze(squeeze(squeeze(reshape(output(:,9),nP,nT,nm))));
Archer.rhow=1e3*squeeze(squeeze(squeeze(reshape(output(:,10),nP,nT,nm))));
Archer.alphas=squeeze(squeeze(squeeze(reshape(output(:,11),nP,nT,nm))));
Archer.alphaw=squeeze(squeeze(squeeze(reshape(output(:,12),nP,nT,nm))));
Archer.betas=squeeze(squeeze(squeeze(reshape(output(:,13),nP,nT,nm))));
Archer.betaw=squeeze(squeeze(squeeze(reshape(output(:,14),nP,nT,nm))));
Archer.Cps=1e3*squeeze(squeeze(squeeze(reshape(output(:,15),nP,nT,nm))));
Archer.Cpw=1e3*squeeze(squeeze(squeeze(reshape(output(:,16),nP,nT,nm))));
Archer.Go=1e3*squeeze(squeeze(squeeze(reshape(output(:,17),nP,nT,nm))));
Archer.Ho=1e3*squeeze(squeeze(squeeze(reshape(output(:,18),nP,nT,nm))));
Archer.So=1e3*squeeze(squeeze(squeeze(reshape(output(:,19),nP,nT,nm))));
Archer.Aphi=squeeze(squeeze(squeeze(reshape(output(:,20),nP,nT,nm))));
Archer.Av=squeeze(squeeze(squeeze(reshape(output(:,21),nP,nT,nm))));
Archer.Ac=squeeze(squeeze(squeeze(reshape(output(:,22),nP,nT,nm))));
Archer.Ah=squeeze(squeeze(squeeze(reshape(output(:,23),nP,nT,nm))));
Archer.diel=squeeze(squeeze(squeeze(reshape(output(:,24),nP,nT,nm))));
Archer.Gw=1e3*squeeze(squeeze(squeeze(reshape(output(:,25),nP,nT,nm))));
Archer.Sw=1e3*squeeze(squeeze(squeeze(reshape(output(:,26),nP,nT,nm))));

% remove singletons from arrays 
if ((Pflg) && (Tflg) && (mflg))
    mm=squeeze(squeeze(squeeze(mm(1,1,1))));
    tm=squeeze(squeeze(squeeze(tm(1,1,1))));
elseif (not(Pflg) && (Tflg) && (mflg))
    mm=squeeze(squeeze(mm(:,1,1)));
    tm=squeeze(squeeze(tm(:,1,1)));
elseif (not(Pflg) && not(Tflg) && (mflg))
    mm=squeeze(squeeze(mm(:,:,1)));
    tm=squeeze(squeeze(tm(:,:,1)));
elseif ((Pflg) && not(Tflg) && not(mflg))
    mm=squeeze(squeeze(mm(1,:,:)));
    tm=squeeze(squeeze(tm(1,:,:)));
elseif (not(Pflg) && not(Tflg) && (mflg))
    mm=squeeze(squeeze(mm(:,:,1)));
    tm=squeeze(squeeze(tm(:,:,1)));
elseif (not(Pflg) && (Tflg) && not(mflg))
    mm=squeeze(squeeze(mm(:,1,:)));
    tm=squeeze(squeeze(tm(:,1,:)));
elseif ((Pflg) && (Tflg) && not(mflg))
    mm=squeeze(squeeze(mm(1,1,:)));
    tm=squeeze(squeeze(tm(1,1,:)));
end

% use definitions to calculate Gibbs energy from the parts
% log of water activity
lnaw=-mm.*Archer.phi*2/omega;

% G per kilogram of water
Archer.G= Archer.Gw + omega*R*tm.*lnaw + mm.*(Archer.Go +  R*2*tm.*log(eps+mm.*Archer.gam));

% also G per kilogram of solution
fac=((1+MW*mm)).^-1;
Archer.Gs=fac.*Archer.G;


