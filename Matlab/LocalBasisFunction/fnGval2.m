function  [rho,vel,G,Cp,alpha,S,U,H,K,Kp,mus,muw]=fnGval(sp,input,M)
% function to return rho,vel, G, Cp, alpha S U H K and Kp for G splines in either (P and T) or (m, P, and T)
%  ALL MKS with P in MPa
%    Usage: [rho,vel,G,Cp,alpha,S,U,H,K,Kp,mu]==G_sp_eval(sp,input)
%          where input is either a npts by (2 or 3) matrix of scatter points in [P  T m] or input is a cell of {P,T,m}  
%                 m in molality, P in MPa and T in K. 
%           rho in kg/m^3, vel in m/s, G in J/Mg Cp in J/kg/K alpha in K^(-1)
%           mu is dG/dm where m is in units determined externally
%  
% JMB 2015  
nw=1000/18.0152;
mu_flg=1;
if nargin==2
    mu_flg=0;
end
%Determine whether gridded or scattered points
% and if this is P-T or P-T-m spline 
if iscell(input) % gridded data
   nd=length(input);  
else % scattered data in columns
   [~,nd]=size(input);
end
if nd==2   % spline in P and T only
    if iscell(input)  % gridded output
        P=input{1};T=input{2};    
        [Pm,Tm]=ndgrid(P,T);
        G=fnval(sp,{P,T});
        d1T=fnval(fnder(sp,[ 0 1]),{P,T});
        d2T=fnval(fnder(sp,[ 0 2]),{P,T});
        d1P=fnval(fnder(sp,[ 1 0]),{P,T});
        dPT=fnval(fnder(sp,[ 1 1]),{P,T});
        d2P=fnval(fnder(sp,[ 2 0]),{P,T});
        d3P=fnval(fnder(sp,[ 3 0]),{P,T});
    else % scatter output
        Tm=input(:,2);
        Pm=input(:,1);
        G=fnval(sp,input')';
        d1T=fnval(fnder(sp,[ 0 1]),input')';
        d2T=fnval(fnder(sp,[0 2]),input')';
        d1P=fnval(fnder(sp,[1 0]),input')';
        dPT=fnval(fnder(sp,[1 1]),input')';
        d2P=fnval(fnder(sp,[2 0]),input')';
        d3P=fnval(fnder(sp,[3 0]),input')';
    end
    mu=[];
else  % spline in P, T and compositions
    if iscell(input) % gridded output
        P=input{1};T=input{2};  m=input{3};  
        [Pm,Tm,mm]=ndgrid(P,T,m);
        G=fnval(sp,{P,T,m});
        d1T=fnval(fnder(sp,[ 0 1 0]),{P,T,m});
        d2T=fnval(fnder(sp,[ 0 2 0]),{P,T,m});
        d1P=fnval(fnder(sp,[ 1 0 0]),{P,T,m});
        dPT=fnval(fnder(sp,[ 1 1 0]),{P,T,m});
        d2P=fnval(fnder(sp,[ 2 0 0]),{P,T,m});
        d3P=fnval(fnder(sp,[ 3 0 0]),{P,T,m});
        dGdm=fnval(fnder(sp,[ 0 0 1]),{P,T,m});
    else % scatter output
        Tm=input(:,2);
        Pm=input(:,1);
        mm=input(:,3);
        G=fnval(sp,input')';
        d1T=fnval(fnder(sp,[ 0 1 0]),input');
        d2T=fnval(fnder(sp,[0 2 0]),input')';
        d1P=fnval(fnder(sp,[1 0 0]),input')';
        dPT=fnval(fnder(sp,[1 1 0]),input')';
        d2P=fnval(fnder(sp,[2 0 0]),input')';
        d3P=fnval(fnder(sp,[3 0 0]),input')';
        dGdm=fnval(fnder(sp,[ 0 0 1]),input')';
    end
end

Cp=-d2T.*Tm;
S=-d1T;
vel=real(sqrt(d1P.^2./(dPT.^2./d2T - d2P))); % MPa-Pa units conversion cancels
rho=1e6*d1P.^(-1);  % 1e6 for MPa to Pa
alpha=1e-6*dPT.*rho; % 1e6 for MPa to Pa
U=G-1e6*Pm./rho+Tm.*S;
H=U-Tm.*S;
K=-d1P./d2P;
Kp=d1P.*d2P.^(-2).*d3P -1;


if mu_flg==1
 f=1+M*mm;
  mus=M*G+f.*dGdm;
  muw=G/nw - 1/nw*f.*mm.*dGdm;
else
    mus=[];
    muw=[];
end