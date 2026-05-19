function  Results=fnG1barval(sp,input,M)
% function to return rho,vel, G, Cp, alpha S U H K and Kp for G splines in either (P and T) or (m, P, and T)
%  ALL MKS with P in MPa
%    Usage: Results=G_sp_eval(sp,input, MW) 
%   Results=Result.rho,vel,G,Cp,alpha,S,U,H,Kt,Kp,Ks,mus,muw
%          where input is either a npts by (2 or 3) matrix of scatter points in [P  T m] or input is a cell of {P,T,m}  
%                 m in molality, P in MPa and T in K.  optional MW is molecular
%                 weight in kg/mol - needed for chemical potentials and
%                 partial molar quantities
%           rho in kg/m^3, vel in m/s, G in J/Mg Cp in J/kg/K alpha in K^(-1)
%           mu is dG/dm where m is in units determined externally
%  
% JMB 2015  
nw=1000/18.0152;  % number of moles of water in a kilogram of water
mu_flg=1;  

%Determine whether this is P-T or P-T-m spline for gridded or scattered points 
if iscell(input) % gridded data
   nd=length(input);  
else % scattered data in columns
   [~,nd]=size(input);
end


 % spline in T and compositions
if iscell(input) % gridded output
        T=input{1};  m=input{2}; 
        m(m==0)=eps; % add eps to zero concentrations
        mflg=0;
        if mu_flg
        if (m(1)~=eps)
          m=[eps;m(:)];  % add in a zero value in order to calculate apparent quantities and remove it later
          mflg=1;
        end
        end
        [Tm,mm]=ndgrid(T,m);
        G=fnval(sp,{T,m});
        d1T=fnval(fnder(sp,[  1 0]),{T,m});
        d2T=fnval(fnder(sp,[  2 0]),{T,m});
        d3Tm=fnval(fnder(sp,[  2 1]),{T,m});
        dGdm=fnval(fnder(sp,[  0 1]),{T,m});
else % scatter output
        Tm=input(:,1);
        mm=input(:,2); 
        m=mm;
        if mu_flg
          mm0=zeros(size(mm))+eps;
          mflg=1;
        else
            mflg=0;
        end
        G=fnval(sp,input')';
        d1T=fnval(fnder(sp,[ 1 0]),input')';
        d2T=fnval(fnder(sp,[2 0]),input')';
        d3Tm=fnval(fnder(sp,[ 2 1]),input')';
        dGdm=fnval(fnder(sp,[  0 1]),input')';
        % calculate zero concentration derivatives to determine apparent
        % Volume and specific heat
        if mu_flg
           d2T0=fnval(fnder(sp,[ 2 0]),[ Tm(:) mm0(:)]')';
           Cp0=-d2T0.*Tm(:);
        end
end


Cp=-d2T.*Tm;
S=-d1T;

U=G+Tm.*S;
H=U-Tm.*S;

if iscell(input) % gridded output
  if mu_flg==1   
     f=1+M*mm;
     mus=M*G + f.*dGdm;
     muw=G/nw - 1/nw*f.*mm.*dGdm;
     Cpm=M*Cp - f.* d3Tm.*Tm;
     
     Cpa=(Cp.*f -repmat(squeeze(Cp(:,1)),1,length(m)))./mm;

%        Cpa=Cpa(:,2:end);
% 
%        Cpm=Cpm(:,2:end);
%        G=G(:,2:end);
%        Cp=Cp(:,2:end);
%        S=S(:,2:end);
%   
%        U=U(:,2:end);
%        H=H(:,2:end);
% 
%        mus=mus(:,2:end);
%        muw=muw(:,2:end);
%        f=f(:,2:end);
    
  else
    mus=[];
    muw=[];
    Cpa=[];
    Va=[];
    f=[];
    Vm=[];
    Cpm=[];
  end
else  % Scattered data
   if mu_flg==1   
      f=1+M*mm;
      mus=M*G + f.*dGdm;
      muw=G/nw - 1/nw*f.*mm.*dGdm;
      Cpm=M*Cp - f.* d3Tm.*Tm;
      Cpa=(Cp.*f - Cp0)./mm;
   else
     mus=[];
     muw=[];
     Cpa=[];
     f=[];
     Cpm=[];
   end
end

Results.Cpa=Cpa;
Results.Cp=Cp;
Results.G=G;
Results.S=S;
Results.U=U;
Results.H=H;
Results.mus=mus;
Results.muw=muw;
Results.f=f;
Results.Cpm=Cpm;