function  Results=fnGval(sp,input)
% function to return properties for G splines in either (P and T) or (m, P, and T)
%  Mostly MKS with P in MPa  
%       exception:  partial molal and apparent volumes in cc/mol 
%
% Usage: Results=fnGval(sp,PTm) 
%          where input is in standard PTm format of a matrix for scattered points
%          or a structure for gridded values. 
%     examples:
%        scattered output with two points of analysis:
%                  PTm=[ 1000 300 1
%                        2000 300 1];
%        gridded output with vectors of PTm:
%                 PTm={0:1000, 300:400, 0:6};
%
% can handle dimensionless T as log(T/Tc) if spline contains field "Tc" 
%
% Note: solute activity coefficient requires a definition of the "standard state". 
%   If the spline field "Go" is set to zero, the standard state is set to zero. If sp.Go is a univarient spline in temperature, 
%  the standard state is calculated and determinations of the activity coefficient and Gex are valid.
%  
% Results = 
%   struct with fields: 
%  Summary of the input:
%         P:   Pressure used (MPa)
%         T:   Temperatures used (K)
%         m:   concentration (kg/mol)
%         x:   concentration as mole fraction 
%         f:   factor to convert from kg of solution to solution containing a kg of water
%
% The following extensive quantities are all per kg of solution
%         G:   Gibbs energy (J/kg) 
%         S:   Entropy (J/K/kg) 
%         U:   Internal Energy (J/kg) 
%         H:   Enthalpy (J/kg) 
%         A:   Helmholz energy (J/kg) 
%         F:   Helmholz energy (J/kg)  
%
% Other standard thermodynamic properties
%       rho:   Density (kg/m^3)
%        Cp:   constant pressure specific heat (J/kg/K)
%        Cv:   constant volume specific heat (J/kg/K)  
%        Kt:   constant T bulk modulus (MPa)
%        Ks:   constant S bulk modulus (MPa)
%        Kp:   pressure derivative of Kt
%     alpha:   thermal expansivity (1/K)
%       vel:   sound speed (m/s)
%        Js:   dT/dPs (K/MPa) experimentally accessible quantity
%     gamma:   Gruneisen parameter
%
% The following are solution properties:
%        Vm:   partial molal volume of solute (cc/mol)
%        Vw:   partial molal volume of solvent (cc/mol)
%        Va:   Apparent volume (cc/mol)
%       Vex:   excess volume (cc/mol)
%       Cpa:   Apparent specific heat (J/K/mol)
%       Cpm:   partial molal specific heat (J/K/mol)
%       mus:   chemical potential of solute (J/mol)
%       muw:   chemical potential of solvent (J/mol)
%       phi:   solvent property
%        aw:   activity of solvent 
%
% The following depend on having included a description of the standard state as a function of temperature at 1 bar    
%       gam:   activity coefficient (valid only if standard state is correctly parameterized)
%       Gex:   excess Gibbs energy (J/mol) (valid only if standard state is correctly parameterized)
%
%   The code uses Pitzer's separation of excess and concentration dependent
%   components. In his parameterization, the part that is not excess is
%   also not ideal mixing. 
%
%      mask:   mask if used (NaN values are interpolated to the PTm points
%      of evaluation.)  This is useful to mask off regions that are clearly
%      not well represented.
%
% The following was added in 2023 to test behavior in solutions of a quantity that ought to be "small" 
%      dCpdPd2m:  second concentration derivative of the pressure derivative of Cp
% 
%  Notes regarding apparent and partial molal properties:  The apparent
%  qualtities are "well behaved"  because they require just an evaluation of the 
%  difference between the solvent and the solvent plus a little solute. 
%  In contrast, partial molal properties
%  involve derivatives of the properties.  As concentration goes to zero
%  one is dividing a derivative by a small number.  The typical first knot
%  in concentration for splines is O(10^-5). Experience shows that the
%  partial properties begin to behave poorly for concentrations smaller
%  than O(10^-4).  In the current code, "infinite dilution" is taken to be
%  2x10-4 in order to avoid the low concentration problems.  This is a
%  topic in need of constant review and reconsideration
%
% JMB 2015-2024

nw=1000/18.01528;  % number of moles of water in a kilogram of water
nw=1000/18.015268; % Tilner Roth and Friend (NIST) value
R=8.3144;  %kJ/mol/K

%Determine whether this is P-T or P-T-m spline for gridded or scattered points 
mu_flg=1;
if length(sp.knots)==2
    mu_flg=0;
end

flgTc=0;
%determine whether dimensionless T is used
if isfield(sp,'Tc')
    flgTc=1;
    Tc=sp.Tc;
end
    
% determine whether scattered or gridded 
if iscell(input) % gridded data
   nd=length(input);  
else % scattered data in columns
   [~,nd]=size(input);
end

%handle P T representations
if nd==1
    error('Check input PTm: fnGval finds one column of input. Input needs to be column centric: [P(:) T(:) optional m(:)]')
elseif nd==2   % spline in P and T only
    if iscell(input)  % gridded output
        P=input{1};T=input{2};    
        nP=length(P);
        nT=length(T);
        [Pm,Tm]=ndgrid(P,T);
        if flgTc
            tau=log(T/Tc);
            %taum=log(Tm/Tc);
            tautm=Tm.^-1;
            tau2tm=-Tm.^-2;   
        else
            tau=T;
        end

        G=sp_val(sp,{P,tau});
        d1T=sp_val(sp,[ 0 1],{P,tau});
        d2T=sp_val(sp,[ 0 2],{P,tau});
        d1P=sp_val(sp,[ 1 0],{P,tau});
        dPT=sp_val(sp,[ 1 1],{P,tau});
        d2P=sp_val(sp,[ 2 0],{P,tau});
        d3P=sp_val(sp,[ 3 0],{P,tau});
    else % scatter output
        Tm=input(:,2);
        Pm=input(:,1);
        if flgTc
            tau=log(Tm/Tc);
            tautm=Tm.^-1;
            tau2tm=-Tm.^-2;   
        else
            tau=Tm;
        end
        input=[input(:,1) tau];

        G=sp_val(sp,input);
        d1T=sp_val(sp,[ 0 1],input);
        d2T=sp_val(sp,[0 2],input);
        d1P=sp_val(sp,[1 0],input);
        dPT=sp_val(sp,[1 1],input);
        d2P=sp_val(sp,[2 0],input);
        d3P=sp_val(sp,[3 0],input);
    end

elseif nd==3  % spline in P, T and concentration  
    if (isfield(sp,'MW'))
        M=sp.MW;
 %       default.  In the  future water may not be  the solvent and the
 %       first enty gives the MW for the solvent.
        if length(M)==2
           M=M(2);
        end
    else
        error('a compositional LBF needs the molecular weight set in the structure')
    end
    if (isfield(sp,'nu'))
        nu=sp.nu;
    else
        error('a compositional LBF needs the number of ions in solution (nu) set in the structure')
    end
    if (isfield(sp,'Go'))
        Goin=sp.Go;        
    else
        Goin=0;
    end
    if iscell(input) % gridded output
        P=input{1};T=input{2};  m=input{3}; 
        x=m/(m+nw);
        m(m==0)=eps; % add eps to zero concentrations to avoid divide by zero
        mflg=0;
        if mu_flg
          nm=length(m);
        if (m(1)>eps)
          m=[2e-4;m(:)];  % add in a near zero value in order to calculate partial quantities and remove it later
          mflg=1;
        end
        end   

        [Pm,Tm,mm]=ndgrid(P,T,m);

        nP=length(P);
        nT=length(T);
        
        if flgTc
            % gotta determine Tc for each x (which is entered as m). The
            % problem here is that the spline evaluations become
            % "scattered" rather than being 'gridded"
            tau=log(T/Tc);
            tautm=Tm.^-1;
            tau2tm=-Tm.^-2;   
        else
            tau=T;
        end

        G=sp_val(sp,{P,tau,m});
        d1T=sp_val(sp,[ 0 1 0],{P,tau,m});
        d1Tm=sp_val(sp,[ 0 1 1],{P,tau,m});
        d2T=sp_val(sp,[ 0 2 0],{P,tau,m});
        d2Tm=sp_val(sp,[ 0 2 1],{P,tau,m});
        dP2T2m=sp_val(sp,[1 2 2],{P,tau,m});
        d1P=sp_val(sp,[ 1 0 0],{P,tau,m});
        dPm=sp_val(sp,[ 1 0 1],{P,tau,m});
        dPT=sp_val(sp,[ 1 1 0],{P,tau,m});
        d2P=sp_val(sp,[ 2 0 0],{P,tau,m});
        d3P=sp_val(sp,[ 3 0 0],{P,tau,m});
        dGdm=sp_val(sp,[ 0 0 1],{P,tau,m});
    else % scattered output
        Tm=input(:,2);
        Pm=input(:,1);
        mm=input(:,3); 
        if flgTc
            tau=log(Tm/Tc);
            tautm=Tm.^-1;
            tau2tm=-Tm.^-2;   
        else
            tau=Tm;
        end
        if(isstruct(Goin))
            Gss=sp_val(Goin,Tm);  % using a univarient spline for the standard state vs T
        else
            Gss=0;
        end
     
        m=mm;
        x=mm/(mm+nw);
        if mu_flg
          mm0=2e-4*ones(size(mm)); % the "fix" to approximate infinite dilution behavior
          mflg=1;
        else
            mflg=0;
        end
        input=[input(:,1) tau input(:,3)];
        G=sp_val(sp,input);
        d1T=sp_val(sp,[ 0 1 0],input);
        d1Tm=sp_val(sp,[ 0 1 1],input);
        d2T=sp_val(sp,[0 2 0],input);
        d2Tm=sp_val(sp,[0 2 1],input);
        dP2T2m=sp_val(sp,[1 2 2],input);
        d1P=sp_val(sp,[1 0 0],input);
        dPm=sp_val(sp,[1 0 1],input);
        dPT=sp_val(sp,[1 1 0],input);
        d2P=sp_val(sp,[2 0 0],input);
        d3P=sp_val(sp,[3 0 0],input);
        dGdm=sp_val(sp,[ 0 0 1],input);
        
%        calculate zero concentration derivatives to determine apparent Volume and specific heat
        if mu_flg
           d2T0=sp_val(sp,[0 2 0],[Pm(:) tau(:) mm0(:)]);
           d1T0=sp_val(sp,[0 1 0],[Pm(:) tau(:) mm0(:)]);
           d1P0=sp_val(sp,[1 0 0],[Pm(:) tau(:) mm0(:)]);
           dPm0=sp_val(sp,[1 0 1],[Pm(:) tau(:) mm0(:)]);
           Gw=sp_val(sp,[Pm(:) tau(:) mm0(:)]);
           if flgTc
              d2T0=d2T0.*tautm.^2+d1T0.*tau2tm;
           end
           V0=1e-6*d1P0;  % 1e6 originates in the MPa to Pa conversion
           Cp0=-d2T0.*Tm(:);
        end
    end
end

if flgTc
    d2T=d2T.*tautm.^2+d1T.*tau2tm;    
    d1T=d1T.*tautm;   
    dPT=dPT.*tautm;
    if mu_flg
     d2Tm=d2Tm.*tautm.^2+d1Tm.*tau2tm;
     d1Tm=d1Tm.*tautm;
    end
end

Cp=-d2T.*Tm;
Cv= Cp + Tm.*dPT.^2./d2P;
S=-d1T;
rho=1e6*d1P.^(-1);  % 1e6 for MPa to Pa

vel=real(sqrt(d1P.^2./(dPT.^2./d2T - d2P))); % MPa-Pa units conversion cancels. 
% Problem with poor parameters that might produce an imaginary sound
% speed is overcome by reporting just the real component.  

Ks=rho.*vel.^2/1e6; % 1e6 for MPa to Pa
alpha=1e-6*dPT.*rho; % 1e6 for MPa to Pa

U=G-1e6*Pm./rho+Tm.*S;% 1e6 for MPa to Pa
A=U-Tm.*S;
H=G+Tm.*S;
Kt=-d1P./d2P;
Kp=d1P.*d2P.^(-2).*d3P -1;
Js=-dPT./d2T;  % an experimentally accessible quantity: the T change with adiabatic pressure change

if iscell(input) % gridded output
  if mu_flg==1   
     f=1+M*mm;
     ntot=nw+mm;
     x=mm./ntot;
     V=rho.^-1;     
     % Calculate chemical potentials
     mus=M*G + f.*dGdm;
     muw=G/nw-f.*mm.*dGdm/nw;
     % deal with standard state:
     if(isstruct(Goin))
       Gss=fnval(Goin,T);  % using a univarient spline in T for the standard state at 1 bar
       % dGss is change in chemical potential at near zero concentration
       % from 1 bar to high P  - this should be a reasonable determination
       % for the standard state cheange at high P
       G2=sp_val(sp,{P,tau,1e-3});
       dGdm2=sp_val(sp,[ 0 0 1],{P,tau,1e-3});% dG/dm for standard state concentration of m=1
       G1b=sp_val(sp,{.1,tau,1e-3});
       dGdm1=sp_val(sp,[ 0 0 1],{.1,tau,1e-3});% dG/dm for standard state concentration of m=1
       dGss=(M*G2+dGdm2)-(M*G1b+dGdm1); 
       Gss=Gss+dGss;
     else
        Gss=0;
        warning('standard state set to zero')
      end

     % calculate other partial  and apparent molal properties
     Vm=(M*d1P +f.*dPm);
     Vw=d1P/nw-f.*mm.*dPm/nw; 
     Cpm=Cp*M - f.* d2Tm.*Tm;
     Cpa=(Cp.*f -repmat(squeeze(Cp(:,:,1)),1,1,length(m)))./mm;
     Va=1e6*(V.*f - repmat(squeeze(V(:,:,1)),1,1,length(m)))./mm;
     Vex= Va-Vm(:,:,1);  
     phi=nw*(G(:,:,1)/nw-muw)./mm./Tm/R/nu;
     aw=exp(-mm.*phi*2/nw);
 % Below depends on having a description of the "standard state"
 % two equivalent versions of the calculations for Gex and gamma are possible
        % Gex= G.*f-G(:,:,1)+mm.*(-Gss+R*nu*Tm.*(1-log(mm)))  
        % gam=exp(Gex./mm./Tm/R/nu-1+phi);
     gam=exp(1/R*1/nu*(mus-Gss)./Tm - log(mm));
     Gex=R*nu*Tm.*mm.*(log(gam)+(1-phi));
     dCpdPd2m=-dP2T2m.*Tm;
     if mflg  % remove the added zero concentration point if added above
       mm=mm(:,:,2:end);
       Pm=Pm(:,:,2:end);
       Tm=Tm(:,:,2:end);
       Va=Va(:,:,2:end);
       Cpa=Cpa(:,:,2:end);
       Vm=Vm(:,:,2:end);
       Vw=Vw(:,:,2:end);
       Cpm=Cpm(:,:,2:end);
       G=G(:,:,2:end);
       rho=rho(:,:,2:end);
       Cp=Cp(:,:,2:end);
       Cv=Cv(:,:,2:end);
       S=S(:,:,2:end);
       vel=vel(:,:,2:end);
       Ks=Ks(:,:,2:end);
       alpha=alpha(:,:,2:end);
       Js=Js(:,:,2:end);
       U=U(:,:,2:end);
       H=H(:,:,2:end);
        A=A(:,:,2:end);
       Kt=Kt(:,:,2:end);
       Kp=Kp(:,:,2:end);
       mus=mus(:,:,2:end);
       muw=muw(:,:,2:end);
       f=f(:,:,2:end);
       gam=gam(:,:,2:end);
       phi=phi(:,:,2:end);
       Vex=Vex(:,:,2:end);
       Gex=Gex(:,:,2:end);
       aw=aw(:,:,2:end);
         x=x(:,:,2:end);
      
     end
  else
    mus=[];
    muw=[];
    Cpa=[];
    Va=[];
    f=[];
    Vm=[];
    Cpm=[];
    gam=[];
    phi=[];
    Gss=[];
    Vex=[];
    Gex=[];
    aw=[];
  end
else  % Scattered data  -
   if mu_flg==1   
      f=1+M*mm;
      V=rho.^-1;
      mus=M*G + f.*dGdm;
      muw=G/nw - 1/nw*f.*mm.*dGdm;
      dCpdPd2m=-dP2T2m.*Tm;
      Vm=M*d1P +f.*dPm;
      Vmo=M*d1P0 +(1+M*eps)*dPm0;
      Vw=d1P/nw-f.*mm.*dPm/nw; 
      Cpm=Cp*M - f.* d2Tm.*Tm;
      Cpa=(Cp.*f - Cp0)./mm;
      Va=1e6*(V.*f - V0)./mm;
      Vex=Va-Vmo;  % need infinite dilution Vm
      phi=nw*(Gw/nw-muw)./mm./Tm/R/nu; % Need G for water
      aw=exp(-mm.*phi*2/nw);
       if(isstruct(Goin))
       Gss=fnval(Goin,Tm);  % using a univarient spline in T for the standard state at 1 bar
       % dGss is change in chemical potential at near zero concentration
       % from 1 bar to high P  - this should be a reasonable determination
       % for the standard state cheange at high P
       G2=sp_val(sp,[Pm(:) tau(:) 1e3*ones(length(mm0(:)),1)]);
       dGdm2=sp_val(sp,[ 0 0 1],[Pm(:) tau(:) 1e3*ones(length(mm0(:)),1)]);% dG/dm for standard state concentration of m=1
       G1b=sp_val(sp,[1e-1*ones(length(mm0(:)),1) tau(:) 1e3*ones(length(mm0(:)),1)]);
       dGdm1=sp_val(sp,[ 0 0 1],[1e-1*ones(length(mm0(:)),1) tau(:) 1e3*ones(length(mm0(:)),1)]);% tau(:) 1e3*ones(length(mm0(:)),1)]);% dG/dm for standard state concentration of m=1
       dGss=(M*G2+dGdm2)-(M*G1b+dGdm1); 
       Gss=Gss+dGss;
     else
        Gss=0;
        warning('standard state set to zero')
       end
      gam=exp(1/R*1/nu*(mus-Gss)./Tm - log(mm));
      Gex=R*nu*Tm.*mm.*(log(gam)+(1-phi));
   else
     mus=[];
     muw=[];
     Cpa=[];
     Va=[];
     f=[];
     Vm=[];
     Cpm=[];
     gam=[];
     phi=[];
     Gss=[];
     Vex=[];
     Gex=[];
     aw=[];
   end
end

Results.G=G;
Results.S=S;
Results.U=U;
Results.H=H;
Results.A=A;
Results.F=A;
Results.rho=rho;
Results.Cp=Cp;

Results.Cv=Cv;
Results.Kt=Kt;
Results.Ks=Ks;
Results.Kp=Kp;
Results.alpha=alpha;
Results.vel=vel;
Results.Js=Js;
Results.gamma=1e6*Results.alpha.*Results.Kt./Results.rho./Results.Cv;
if iscell(input)
   Results.P=P;
    Results.T=T; 
else
    Results.P=Pm;
    Results.T=Tm;
end
if mu_flg
    Results.Va=Va;
    Results.Cpa=Cpa;
    Results.dCpdPd2m=dCpdPd2m;
    Results.mus=mus;
    Results.muw=muw;
    Results.f=f;
    Results.m=mm;
    Results.x=x;
    Results.Vm=Vm;
    Results.Vw=Vw;
    Results.Cpm=Cpm;
    Results.gam=gam;
    Results.phi=phi;
    Results.Vex=Vex;
    Results.Gex=Gex;
    Results.aw=aw;    
end

if isfield(sp,'mask')  % create output mask for the input independent variables - NaN for anything not in the defined domain
    mask=sp.mask;
    d=length(size(mask));
    if d==2
        pg=sp.PTm_mask{1};
        tg=sp.PTm_mask{2};
        mask=mask'; 
        if iscell(input) % gridded output
           Results.mask=reshape(interp2(pg,tg,mask,Pm(:),Tm(:)),nP,nT);   
        else
           Results.mask=interp2(pg,tg,mask,Pm(:),Tm(:));   
        end       
    elseif d==3
        pg=sp.PTm_mask{1};
        tg=sp.PTm_mask{2};
        mg=sp.PTm_mask{3};
        mask=permute(mask,[2 1 3]);
        if iscell(input) % gridded output
           Results.mask=reshape(interp3(pg,tg,mg,mask,Pm(:),Tm(:),mm(:)),nP,nT,nm);   
        else
           Results.mask=interp3(pg,tg,mg,mask,Pm(:),Tm(:),mm(:));   
        end
    end
end

function out=sp_val(sp,derv,x)
% function sp_val allows use of LBF representations without
% installation of the MATLAB optional Curevefitting Toolbox .  
% sp_val substitutes for the combinations of fnval and fnder plus various utility 
% functions in the spline toolbox.  sp_val is "self contained" and calls no other
% specialized functions. The algorithms are those given by de Boor.  
%
%  Usage:   out = sp_val(LBF,derv,x) or = sp_val(LBF,x) if derivatives are not needed
%
%  where:    "out" is a vector or matrix of values (depending on the form of"x"
%            "LBF" is a b spline (tensor or univarient)
%            "derv" is vector of requested derivatives (one element per degree of freedom)
%            "x" is either a cell containing the independent variables (gridded output)
%                 or an array of scattered data.
%             sp_val(LBF,[0 0 ..],x) is equivalent to sp_val(LBF,x)
%  example: for an LBF in pressure (0-100 MPa), temperature (250-350 K) and concentration (0-1 M),
%           PTm={0:10:100,250:10:350,0:.1:1} and for first derivatives of
%           all independent variables:
%                  out=sp_val(LBF,[1 1 1],PTm);
%
% no error checking is included - user must understand how to use this function. 
%
% a major limitation of this function is that scattered data are calculated
% in a loop using the b-form splines (parfor speeds it up).  For scattered data, 
% the Mathworks provided fn2fm convert the b-splines to a pp-form and then does the 
% calculation using the pp-form. To do this requires more implementation of various
% spline manipulations. The slower method used here is about 3-5 times slower.
% 7/2019


if nargin==2,x=derv;derv=zeros(size(x(1,:)));end

   sp=fndr(sp,derv);  
   out = spvl(sp,x);
end

end
       
function out = spvl(sp,x)
% based on SPVAL that evaluates a function in B-form.

if iscell(sp.knots)  % we are dealing with a tensor product spline
   [t,a,n,~,d] = spbk(sp); m = length(t);
   %nd=length(x);
   if iscell(x)  % evaluation on a mesh
     v = a; sizev = [d,n]; nsizev = zeros(1,m);
      for i=m:-1:1
         nsizev(i) = length(x{i}(:));
         v = reshape(...
         spv(spmk(t{i},reshape(v,prod(sizev(1:m)),sizev(m+1))), ...
                 x{i}),   [sizev(1:m),nsizev(i)]);
         sizev(m+1) = nsizev(i);
         if m>1
            v = permute(v,[1,m+1,2:m]); sizev(2:m+1) = sizev([m+1,2:m]);
         end
      end
      if d>1
         out = reshape(v,[d,nsizev]);
      else
         out = reshape(v,nsizev);
      end
 
   else          % evaluation at scattered points;
                 % this will eventually be done directly here.
%      if iscell(spline.knots)   % we are dealing with a multivariate spline
% 
%    [t,a,n,k,d] = spbrk(spline);
%    m = length(k);
%    coefs = a; sizec = [prod(d),n]; % size(coefs);
%    for i=m:-1:1
%       ppi = sp2pp1(spmak(t{i},reshape(coefs,prod(sizec(1:m)),n(i))));
%       breaks{i} = ppi.breaks;  sizec(m+1) = ppi.pieces*k(i);
%       coefs = reshape(ppi.coefs,sizec);
%       if m>1
%          coefs = permute(coefs,[1,m+1,2:m]); sizec = sizec([1,m+1,2:m]);
%       end
%    end
%    pp = ppmak(breaks,coefs,sizec);
%       
% else
%    pp = sp2pp1(spline);
% end
            
      [nd,~]=size(x);          
      out=zeros(nd,1);
      
      parfor ii=1:nd
         sizev = [d,n]; nsizev = ones(1,m);
         v = a;
      for i=m:-1:1
 
         v = reshape(...
         spv(spmk(t{i},reshape(v,prod(sizev(1:m)),sizev(m+1))), ...
                 x(ii,i)),   [sizev(1:m),1]);
         sizev(m+1) = 1;
         if m>1
            v = permute(v,[1,m+1,2:m]); sizev(2:m+1) = sizev([m+1,2:m]);
         end
      end
      if d>1
         v = reshape(v,[d,nsizev]);
      else
         v = reshape(v,nsizev);
      end
      out(ii)=v;
      end
                                   
   end
else                 % we are dealing with a univariate spline 

   v = spv(sp,x);
   out=v;
end
end

function v = spv(sp,x)
% based on SPVAL1 that evaluate a univariate function in B-form.

[mx,nx] = size(x); lx = mx*nx; xs = reshape(x,1,lx);
%  Take apart spline:
[t,a,n,k,d] = spbk(sp);
if lx==0, v = zeros(d,0); return, end

%  Otherwise, augment the knot sequence so that first and last knot each
%  has multiplicity  >= K . (AUGKNT would not be suitable for this since
%  any change in T must be accompanied by a corresponding change in A.)

index = find(diff(t)>0); addl = k-index(1); addr = index(end)-n;
if ( addl>0 || addr>0 )
   npk = n+k; t = t([ones(1,addl) 1:npk npk(ones(1,addr))]);
   a = [zeros(d,addl) a zeros(d,addr)];
   n = n+addl+addr;
end

% For each data point, compute its knot interval:
   [~, index] = histc(-xs,[-inf,-fliplr(t(k+1:n)),inf]);
   NaNx = find(index==0); index = max(n+1-index,k);
if ~isempty(NaNx), index(NaNx) = k; end

% Now, all is ready for the evaluation.
if  k>1  % carry out in lockstep the first spline evaluation algorithm
         % (this requires the following initialization):
   dindex = reshape(repmat(index,d,1),d*lx,1);
   tx =reshape(t(repmat(2-k:k-1,d*lx,1)+repmat(dindex,1,2*(k-1))),d*lx,2*(k-1));
   tx = tx - repmat(reshape(repmat(xs,d,1),d*lx,1),1,2*(k-1));
   dindex = reshape(repmat(d*index,d,1)+repmat((1-d:0).',1,lx),d*lx,1);
   b = repmat(d*(1-k):d:0,d*lx,1)+repmat(dindex,1,k);
   a = a(:); b(:) = a(b);

   % (the following loop is taken from deBoor SPRPP)

   for r = 1:k-1
      for i = 1:k-r
         b(:,i) = (tx(:,i+k-1).*b(:,i)-tx(:,i+r-1).*b(:,i+1)) ./ ...
                  (tx(:,i+k-1)    -    tx(:,i+r-1));
      end
   end
   v = reshape(b(:,1),d,lx);
else     % the spline is piecewise constant, hence ...
   v = a(:,index);
   if ~isempty(NaNx), v(:,NaNx) = NaN; end
end

% Finally, zero out all values for points outside the basic interval:
index = find(x<t(1)|x>t(n+k));
if ~isempty(index)
   v(:,index) = zeros(d,length(index));
end
v = reshape(v,d*mx,nx);
end
  
   
function spline = spmk(knots,coefs,sizec)
%based on deBoor SPMAK that puts together a spline in B-form.

if nargin<3
 sizec = size(coefs);
end
m = 1;
if iscell(knots)
    m = length(knots);
end

if length(sizec)==m  % coefficients of a scalar-valued function
    sizec = [1 sizec];
end

% convert ND-valued coefficients into vector-valued ones, retaining the
% original size in SIZEVAL, to be temporaily stored in SP.DIM .
sizeval = sizec(1:end-m);
sizec = [prod(sizeval), sizec(end-m+(1:m))];
coefs = reshape(coefs, sizec);

if iscell(knots) % we are putting together a tensor-product spline
    [knots,coefs,k,sizec] = chkt(knots,coefs,sizec);
else            % we are putting together a univariate spline
    [knots,coefs,k,sizec] = chkt({knots},coefs,sizec);
    knots = knots{1};
end

spline.form = 'B-';
spline.knots = knots;
spline.coefs = coefs;
spline.number = sizec(2:end);
spline.order = k;
spline.dim = sizeval;
end

function [knots,coefs,k,sizec] = chkt(knots,coefs,sizec)
%based on CHCKKNT: checks knots and omits trivial B-splines
for j=1:length(sizec)-1
    n = sizec(j+1);
    k(j) = length(knots{j})-n;

       % make sure knot sequence is a row matrix:
    knots{j} = reshape(knots{j},1,n+k(j));
    % throw out trivial B-splines:
    index = find(knots{j}(k(j)+(1:n))-knots{j}(1:n)>0);
    if length(index)<n
        oldn = n;
        n = length(index);
        knots{j} = reshape(knots{j}([index oldn+(1:k(j))]),1,n+k(j));
        coefs = ...
            reshape(coefs, [prod(sizec(1:j)),sizec(j+1),prod(sizec(j+2:end))]);
        sizec(j+1) = n;
        coefs = reshape(coefs(:,index,:),sizec);
    end
end
end


function fprime = fndr(f,dorder)
%based on FNDER:  differentiates a function by differencing the coefficients (see de Boor).
sizeval = f.dim;
if length(sizeval)>1, f.dim = prod(sizeval); end
if nargin<2, dorder=1; end

   [knots,coefs,n,~,d]=spbk(f);
   if iscell(knots)       % the function is multivariate
      m = length(knots);
      sizec = [d,n];% size(coefs);
      for i=m:-1:1
         dsp = fndb(spmk(knots{i},...
            reshape(coefs,prod(sizec(1:m)),sizec(m+1))),dorder(i));
         knots{i} = dsp.knots; sizec(m+1) = dsp.number;
         coefs = reshape(dsp.coefs,sizec); 
         if m>1
            coefs = permute(coefs,[1,m+1,2:m]);
            sizec(2:m+1) = sizec([m+1,2:m]);
         end
      end
      fprime = spmk(knots,coefs,sizec);
   else
      fprime = fndb(f,dorder);
   end
if length(sizeval)>1, fprime.dim = sizeval; end
end


function fprime = fndb(f,dorder)
%FNDERB Differentiate a univariate function in B-form.
[t,a,n,k,d]=spbk(f);
if k<=dorder
   fprime=spmk(t,zeros(d,n));
elseif dorder<0    % we are to integrate
   error('integration not implemented')
else
   knew=k-dorder;
   for j=k-1:-1:knew
       % here it is: difference the knots and coefficients
      tt=t(j+1+[0:n])-t(1:n+1); z=find(tt>0); nn=length(z);     
      temp=(diff([zeros(1,d);a.'; zeros(1,d)])).';
      a=temp(:,z)./repmat(tt(z)/j,d,1);
      t=[t(z) t(n+2:n+j+1)]; n=nn;
   end
   fprime=spmk(t,a);
end
end


function varargout = spbk(sp)
    varargout = {sp.knots,sp.coefs, sp.number, sp.order, sp.dim};
end


function pp = s2p1(spline)
%  Take apart the  spline

[t,a,n,k,d] = spbk(spline);

%  and augment the knot sequence so that first and last knot each have
%  multiplicity  k .

index = find(diff(t)>0); addl = k-index(1); addr = index(end)-n;
if (addl>0||addr>0)
   t = [repmat(t(1),1,addl) t(:).' repmat(t(n+k),1,addr)];
   a = [zeros(d,addl) a zeros(d,addr)];
end

%  From this, generate the pp description.

inter = find( diff(t)>0 ); l = length(inter);
if k>1
   temp = repmat(inter,d,1); dinter = temp(:);
   tx = repmat(2-k:k-1,d*l,1)+repmat(dinter,1,2*(k-1)); tx(:) = t(tx);
   tx = tx-repmat(t(dinter).',1,2*(k-1)); a = a(:);
   temp = repmat(d*inter,d,1)+repmat((1-d:0).',1,l); dinter(:) = temp(:);
   b = repmat(d*(1-k:0),d*l,1)+repmat(dinter,1,k); b(:) = a(b);
   c = srp(tx,b);
else temp = a(:,inter); c = temp(:);
end

%   put together the  pp

pp = ppmk([t(inter) t(inter(end)+1)],c,d);
end


function [v,b] = srp(tx,a)
%SRP Right Taylor coefficients from local B-coefficients.
%
%   [V,B] = SRP(TX,A)
%
%   uses knot insertion to derive from the B-spline coefficients
%   A(.,:) relevant for the interval  [TX(.,k-1) .. TX(.,k)]  (with
%   respect to the knot sequence  TX(.,1:2k-2) )  the polynomial
%   coefficients V(.,1:k) relevant for the interval  [0 .. TX(.,k)] .
%   Here,    [ ,k] := size(A) .
%   Also, it is assumed that  TX(.,k-1) <= 0 < TX(.,k) .
%
%   In the process, uses repeated insertion of  0  to derive, in
%   B(.,1:k) , the B-spline coefficients relevant for the interval
%   [0 .. TX(.,k)]  (with respect to the knot sequence
%   [0,...,0,TX(.,k:2*(k-1))]) .
%
%   See also SPLPP.

%   Carl de Boor 25 feb 89
%   Copyright 1987-2008 The MathWorks, Inc. 


k = length(a(1,:)); km1 = k-1; b = a;
for r=1:km1
   for i=1:k-r
      b(:,i) =(tx(:,i+km1).*b(:,i)-tx(:,i+r-1).*b(:,i+1))./...
               (tx(:,i+km1)-tx(:,i+r-1));
   end
end

%  Use differentiation at  0  to generate the derivatives

v = b;
for r=2:k
   factor = (k-r+1)/(r-1);
   for i=k:-1:r
      v(:,i) = (v(:,i) - v(:,i-1))*factor./tx(:,i+k-r);
   end
end

v = v(:,k:-1:1);
end

