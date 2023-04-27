function output=ArcherRardMgSO4EOS(input,H2Oflg)
%this returns the MgSO4 EOS based on Archer and Rard 
% J. Chem. Eng. Data 1998, 43, 791-806
% Isopiestic Investigation of the Osmotic and Activity Coefficients of
% Aqueous MgSO4 and the Solubility of MgSO4â7H2O(cr) at 298.15 K:
% Thermodynamic Properties of the MgSO4 + H2O System to 440 K
%  Usage:  output=ArcherMgSO4EOS(input,H2Oflg)
%        where the output structure contains standard EOS quantities and   activities
%        including Gs (J/kg) Gibbs energy per kg of solution and G per kg of water
%       input:   {T,m} or [T(:) m(:)] with T in (K) and concentration in kg/mol pflg=1 for 200 bars and 0 for 1 bar
%        H2Oflg 'Hill' or 'IAPWS'   Archer and Wang dielectric constant +
%        Hill Water
%        the dielectric constant the second uses Fernandez et al 1997
%  Values still need to be checked against an absolute G based on solubility
% JMB 2017
  
R=8.3144621;
MW=0.120366;
omega=1000/18.0152;
Tr=298.15;
To=1;

% organize input
    if(iscell(input))
        grd_flg=1;
        T=input{1};
        m=input{2};
        m(m==0)=eps;
        %[tm,mm]=ndgrid(T,m);
        nT=length(T);
        nm=length(m); 
    else
        grd_flg=0;
        tm=input(:,1);
        T=tm;
        m=input(:,2)';
        m(m==0)=eps;
        mm=m;
        nT=length(T);
        nm=1;
    end

 % water EOS quantities as column of nT values
 switch H2Oflg
    case 'Hill'
        results = RunArcherNaClFORTRAN({.1,T,0});
        Aphi=results.Aphi(:);
        Cpw=results.Cpw(:);
        muo_w=results.Gw(:); % in J per kg of water
        Av=results.Av(:);  
        Ah=R*T(:).*results.Ah(:);
        Ac=R*results.Ac(:);

    case 'IAPWS'
    
        results = Fernandez1997dielectric({.1,T});
        Aphi=results.Aphi(:);
        Av=results.Av(:);
        Ah=results.Ah(:);
        Ac=results.Ac(:);
        Cpw=results.iapws.Cp(:);
        muo_w=results.iapws.G(:);
        Cpw=results.iapws.Cp(:);
 end

% some constants:
b=1.2;
a1=1.4;
k=30.65;
a2=k*Aphi;
a3=1;


% all compositional quantities as row of nm quantities
fac=(1+MW*m);  % conversion factor from per kilo of solution to per kilo of water
I=4*m;  
sqrtI=sqrt(I);
x1=a1*sqrtI;
x2=a2.*sqrtI;
x3=a3*sqrtI;
 
% The compositional basis functions in rows
g1=2*(1-(1+x1).*exp(-x1))./x1.^2;
g2=2*(1-(1+x2).*exp(-x2))./x2.^2;

g3=4*(6-(6+6*x3+3*x3.^2+x3.^3).*exp(-x3))./x3.^4;
g4=4*(6-(6+6*x3+3*x3.^2+x3.^3-0.5*x3.^4).*exp(-x3))./x3.^4;

h1=2*(1-(1+x1-0.5*x1.^2).*exp(-x1))./x1.^2;
h2=2*(1-(1+x2-0.5*x2.^2).*exp(-x2))./x2.^2;

j1=(1-(1+x1).*exp(-x1))./x1.^2;
j2=(1-(1+x2).*exp(-x2))./x2.^2;
jp2 = -2*x2.^-3.*(1-(1+x2+.5*x2.^2).*exp(-x2));
j2p2 = 6*x2.^-4.*(1-(1+x2+.5*x2.^2+1/6*x2.^3).*exp(-x2));
 
 
% Table 3. Least-Squares Estimated Parameters for the Ion-interaction Model of the Thermodynamic Properties of
% MgSO4(aq) parameter value parameter value parameter valuea
%        beta0                    beta1                       beta2                   C0                       C1
parm=[
    -5.26309458086110e-01     4.79138046549794e+00     1.34733864722672e+03     1.61656666800298e-02     4.11919661612742e-01
     7.98429374952591e-01    -1.95415726851180e+00     6.56287453214204e+02    -4.49163783095252e-02    -1.24114711530086e+00
    -8.62343583047259e+00     2.31818236890809e+01     6.08536936104363e+02     3.31353349416800e-01     1.11400433261952e+01
                        0    -7.52047721452442e-01                        0                        0    -5.55913285968194e-02
     1.91672710497768e+01                         0   -7.36264533444156e+04                        0                        0
                        0                         0    1.61319854373951e+02                        0                        0
    -9.56771923897770e-02                         0    7.31155393709027e+01     1.29599608354122e-02                        0
    ];

% the temperature basis functions  created for viewing ease as rows of T and transposed
T=T(:)';
basis=[ones(1,nT) 
       1e-2*(T-Tr)/To
       1e-5*((T-Tr)/To).^2       
       To*1e2*(T-225).^-1
       1e1*To*(680-T).^-1
       1e3*To*T.^-1      
       1e-7*((T-Tr)/To).^3
       ]';
% first derivative wrt T
basisT=[zeros(1,nT) 
       1e-2/To*ones(1,nT) 
       2e-5*(T-Tr)/To        
       -To*1e2*(T-225).^-2
       1e1*To*(680-T).^-2
       -1e3*To*T.^-2       
       3e-7*((T-Tr)/To).^2
       ]'; 
%second derivative wrt T
basis2T=[
       zeros(1,nT) 
       zeros(1,nT) 
       2e-5/To*ones(1,nT)   
       To*2e2*(T-225).^-3
       2e1*To*(680-T).^-3
       2e3*To*T.^-3  
       6e-7*((T-Tr)/To) 
       ]';
T=T(:);   % need T as column

Cpo=ArcherCp(T);

% multiply basis functions by parameters to get the beta's and c's
BC=basis*parm;
BCT=basisT*parm;
BC2T=basis2T*parm;

% pull out the betas and c's
beta0=BC(:,1);
beta1=BC(:,2);
beta2=BC(:,3);
c0=BC(:,4);
c1=BC(:,5);

% following not yet used
beta0T=BCT(:,1);
beta1T=BCT(:,2);
beta2T=BCT(:,3);
c0T=BCT(:,4);
c1T=BCT(:,5);

beta0TT=BC2T(:,1);
beta1TT=BC2T(:,2);
beta2TT=BC2T(:,3);
c0TT=BC2T(:,4);
c1TT=BC2T(:,5);

% create combinations: 
Bmx=beta0 + beta1.*g1 +beta2.*g2;
Cmx= c0 + c1.*g3;

% for specific heat  - not yet ready to vet this
        Bmx_Cp=beta0TT +2*T.^-1.*beta0T + 2*(beta1TT +2*T.^-1.*beta1T)*j1 + 2*(beta2TT +2*T.^-1.*beta2T).*j2 ...
            + k/R*T.^-2.*sqrtI.*(Ah.*beta2T+1/2*Ac.*beta2).*jp2 ...
              + k^2/R^2/8*T.^-4.*Ah.^2.*beta2.*I.*j2p2;

        Cmx_Cp=(c0TT+2*T.^-1 .*c0T) ...
            + 4*(c1TT+2*T.^-1.*c1T)*((6-(6+6*x3+3*x3.^2+x3.^3).*exp(-x3))./x3.^4);

switch grd_flg
    case 1
        output.Gex=-4*R/b*T.*Aphi*(I.*log(1+b*sqrtI))...
                      +  2*R*T.*Bmx.*m.^2  ...
                           + 4*R*T.*Cmx.*m.^3;
        output.GexDH=4*R/b*T.*Aphi*(I.*log(1+b*sqrtI));
        output.Gexvirial= +  2*R*T.*Bmx.*m.^2  ...
                           + 4*R*T.*Cmx.*m.^3;

        output.phi= 1-4*Aphi*(sqrtI./(1+b*sqrtI))  ...
                        + m.*(beta0 + beta1.*exp(-x1) + beta2.*exp(-x2)) ...
                               + 4*m.^2.*(c0 + c1.*exp(-x3));
         output.phiDH=1-4*Aphi*(sqrtI./(1+b*sqrtI)) ;
         output.phivirial=+ m.*(beta0 + beta1.*exp(-x1) + beta2.*exp(-x2)) ...
                               + 4*m.^2.*(c0 + c1.*exp(-x3));

        output.gamma=exp(-4*Aphi*(sqrtI./(1+b*sqrtI) + 2/b*log(1+b*sqrtI)) ...
            + m.*(2*beta0 + beta1.*h1 +beta2.*h2) ...
                + 2*m.^2.*(3*c0 + c1.*g4));

        output.Cpphi=Cpo+8*Ac.*log(1+b*sqrtI)/2/b ...
                   -2*R*T.^2.*(m.*Bmx_Cp + 2*m.^2.*Cmx_Cp);   
        output.CpDH=8*Ac.*log(1+b*sqrtI)/2/b;
        output.Cpvirial=-2*R*T.^2.*(m.*Bmx_Cp + 2*m.^2.*Cmx_Cp); 
        
        % here is the "non-excess" compositional dependent part of G
        Gnon=R*2*(T.*m).*((log(m)) -1);  %+log(2)/3    
    case 0   

        output.Gex=diag(-4*R/b*T.*Aphi.*(I'.*log(1+b*sqrtI'))...
                      +  2*R*T'.*Bmx.*m.^2  ...
                           + 4*R*T'.*Cmx.*m.^3);


        output.phi= diag(1-4*Aphi.*(sqrtI./(1+b*sqrtI))  ...
                        + m.*(beta0 + beta1.*exp(-x1) + beta2.*exp(-x2)) ...
                               + 4*m.^2.*(c0 + c1.*exp(-x3)));

        output.gamma=diag(exp(-4*Aphi.*(sqrtI./(1+b*sqrtI) + 2/b*log(1+b*sqrtI)) ...
            + m.*(2*beta0 + beta1.*h1 +beta2.*h2) ...
                + 2*m.^2.*(3*c0 + c1.*g4)));

        output.Cpphi=diag(Cpo+8*Ac.*log(1+b*sqrtI)/2/b ...
                           -2*R*T'.^2.*(m.*Bmx_Cp + 2*m.^2.*Cmx_Cp));
        % here is the "non-excess" compositional dependent part of G
        Gnon=diag(R*2*(T'.*m).*((log(m)) -1));  %+log(2)/3
end


                            
output.Cp=Cpw+m.*output.Cpphi;
output.Cps=fac.^-1.*output.Cp;


% 
% G(MgSO47H2O) 10.616 ( 0.066  kJ/mol
% S(MgSO47H2O) 5.24 (1.05  J/K/mol
% G(MgSO46H2O) 8.970(0.073)
% S(MgSO46H2O) -44.15 (0.88)



% need standard state properties of salt and water

% Here is calculatiom for salt based on Pabalan's integration constants

% these values based on total (crystal+solution).
muoos=1e4; % energy per mole at 1 bar and 298.15 K
Sos=5.24;  % delta entropy at 1 bar and 298.15 K for MgSO4*7H2O

% next integrate standard state Cp and Cp/T to get energy at other
% temperatures
Q1=zeros(nT,1);
Q2=Q1;
f2 = @(x) (ArcherCp(x)./x) ; % ACp/T
f1= @(x) ArcherCp(x);  % ACp
for i=1:nT  
        Q1(i)=integral(f1,298.15,T(i));  % integrate from 298 to T
        Q2(i) = integral(f2,298.15,T(i));
end

muo_s=(muoos-Sos*(T-298.15)+Q1-T.*Q2); % this takes a value at 298 and 1 bar and 
mu_s=muo_s + R*2*T.*log(output.gamma.*m);  % arbitrary offsets to slide chemical potential to near zero

lnaw=-m.*output.phi*2/omega;
mu_w=muo_w/omega+R*T.*lnaw;
% this alternate definition has been tested - it works
%G2=omega*mu_w + mm.*mu_s;  % this needs to be equal to the G defined below

% now add them up
Go= muo_w + m.*muo_s;
output.G=output.Gex+Go+Gnon;  % per kg of water
output.Gs=output.G./fac;   % per kg of solution
output.mu_w=mu_w;
output.mu_s=mu_s;
output.Cpw=Cpw;
output.Bmx=Bmx;
output.Cmx=Cmx;


if grd_flg==0
    output.Cp=diag(output.Cp);
    output.Cps=diag(output.Cps);
        output.G=diag(output.G);
            output.Gs=diag(output.Gs);
                output.mu_w=diag(output.mu_w);
                    output.mu_s=diag(output.mu_s);
end


function Cpo=ArcherCp(T)
% input T in K  
Cpo=-295.3 -18.52779*(T-298.15) + 0.0728295*(T.^2-298.15^2)-8.79539e-5*(T.^3-298.15.^3);




