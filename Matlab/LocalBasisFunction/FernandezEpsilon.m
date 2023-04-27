function results=FernandezEpsilon(input)
% A Formulation for the Static Permittivity of Water and Steam at Temperatures from 238
% K to 873 K at Pressures up to 1200 MPa, Including Derivatives and Debye?Hückel
% Coefficients
% D. P. Fernández, A. R. H. Goodwin, Eric W. Lemmon, J. M. H. Levelt Sengers, and R. C. Williams
% Citation: Journal of Physical and Chemical Reference Data 26, 1125 (1997); doi: 10.1063/1.555997
% Usage:results=FernandezEpsilonfunc(input) where 
%                "results" is structure with both dielectric properties and Debye Huckle terms. 
%                "input" is P and T in structure {P,T} or matrix of scattered points [P(:) T(:)] in units of MPa and K
%  IAPWS95 is used for the EOS of water. output has been checked against
%  Tables 12 and 17.  Agreement is better than tens of parts per million for all quantities
%  JMB 2017


% organize input
    if(iscell(input))
        P=input{1};
        T=input{2};
        [pm,tm]=ndgrid(P,T);
        pm=pm(:);
        tm=tm(:);
        nP=length(P);
        nT=length(T);
        npts=length(pm);
    else
        pm=input(:,1);
        tm=input(:,2);
        nP=length(pm);
        nT=1;
        npts=nP;
    end

%Permittivity of free space,  C^2/J/m
eo=(4e-7*pi*299792458^2)^-1;
%Elementary charge, 
e =1.60217733e-19;
%Boltzmann?s constant, J/K
k= 1.380658e-23;
%Avogadro?s number, 
Na= 6.0221367e23;
R=k*Na;
%Molar mass of water, kg/mol
Mw =0.018015268;
%Mean molecular polarizability of water, C^2/J/m^2
a= 1.636e-40;
%Dipole moment of water, C m
u=6.138e-30;
rhoc=322;
Tc=647.096;
dT=1e-7;

parms=[
1     0.978224486826       1    0.25
2    -0.957771379375       1    1
3     0.237511794148       1    2.5
4     0.714692244396       2    1.5
5    -0.298217036956       3    1.5
6    -0.108863472196       3    2.5
7     0.949327488264e-1    4    2
8    -0.980469816509e-2    5    2
9     0.165167634970e-4    6    5
10    0.937359795772e-4    7    0.5
11   -0.123179218720e-9   10   10];
p12 =[0.196096504426e-2   -1.2];

iapws=IAPWS95([pm tm],'P');
iapws2=IAPWS95([pm tm+dT],'P');
rho=iapws.rho;
Kt=iapws.Kt;
alpha=iapws.alpha;
rho2=iapws2.rho;
%Kt2=iapws2.Kt;
alpha2=iapws2.alpha;

tm2=tm+dT;
x=rho/rhoc;
x2=rho2/rhoc;
xp=1/rhoc;
y=Tc*tm.^-1;
y2=Tc*tm2.^-1;
yp=-Tc*tm.^-2;
yp2=-Tc*tm2.^-2;
y12=(tm/228 -1);
y12_2=(tm2/228 -1);
y12p=1/228;

    g=ones(npts,1);
    g2=g;
    dgdt=g;
    dgdt2=g;
    dgdr=g;
    dgdr2=g;
    for i=1:npts
        g(i)=1+sum(parms(:,2).*x(i).^parms(:,3).*y(i).^parms(:,4)) + p12(1)*x(i)*y12(i)^p12(2);
        dgdt(i)=sum(parms(:,2).*x(i).^parms(:,3).*parms(:,4).*y(i).^(parms(:,4)-1).*yp(i)) + y12p*p12(2)*p12(1)*x(i)*y12(i)^(p12(2)-1);
        dgdr(i)=sum(xp*parms(:,3).*parms(:,2).*x(i).^(parms(:,3)-1).*y(i).^parms(:,4)) + p12(1)*xp*y12(i)^p12(2);
        g2(i)=1+sum(parms(:,2).*x2(i).^parms(:,3).*y2(i).^parms(:,4)) + p12(1)*x2(i)*y12_2(i)^p12(2);      
        dgdt2(i)=sum(parms(:,2).*x2(i).^parms(:,3).*parms(:,4).*y2(i).^(parms(:,4)-1).*yp2(i)) + y12p*p12(2)*p12(1)*x2(i)*y12_2(i)^(p12(2)-1);
        dgdr2(i)=sum(xp*parms(:,3).*parms(:,2).*x2(i).^(parms(:,3)-1).*y2(i).^parms(:,4)) + xp*p12(1)*y12_2(i)^p12(2);
    end

A=Na*u^2/eo/k/Mw*rho.*g./tm;
B=Na*a/3/eo/Mw*rho;
A1=A./rho +(A./g).*dgdr;
B1=B./rho;
A2=-A./tm+(A./g).*dgdt;
C=9+2*A+18*B+A.^2+10*A.*B+9*B.^2;
epsilon=(1 + A + 5*B + sqrt(9 + 2*A + 18*B + A.^2 + 10*A.*B + 9*B.^2))./(4-4*B);
dedr_T=((4-4*B).^-1).*(4*B1.*epsilon + A1+5*B1 + 0.5*C.^(-0.5).*(2*A1 +18*B1 + 2*A.*A1 + 10*(A1.*B +A.*B1) + 18*B.*B1));
dedT_rho=((4-4*B).^-1).*(A2 +.5*C.^(-0.5).*A2.*(2+2*A+10*B));
dedP_T=dedr_T.*(rho./Kt);
dedT_P=dedT_rho-dedr_T.*alpha.*rho;

A_phi=1/3*(2*pi*Na*rho).^.5.*(e.^2*(4*pi*eo*k*epsilon.*tm).^-1).^(3/2);
Av=2*R*tm.*A_phi.*(3*dedP_T./epsilon - Kt.^-1);
Ah=-6*R*tm.*A_phi.*(1+tm.*dedT_P./epsilon +tm.*alpha/3);

A=Na*u^2/eo/k/Mw*rho2.*g2./tm2;
B=Na*a/3/eo/Mw*rho2;
A1=A./rho2 +(A./g2).*dgdr2;
B1=B./rho2;
A2=-A./tm2+(A./g2).*dgdt2;
C=9+2*A+18*B+A.^2+10*A.*B+9*B.^2;
epsilon2=(1 + A + 5*B + sqrt(9 + 2*A + 18*B + A.^2 + 10*A.*B + 9*B.^2))./(4-4*B);
dedr_T=((4-4*B).^-1).*(4*B1.*epsilon2 + A1+5*B1 + 0.5*C.^(-0.5).*(2*A1 +18*B1 + 2*A.*A1 + 10*(A1.*B +A.*B1) + 18*B.*B1));
dedT_rho=((4-4*B).^-1).*(A2 +.5*C.^(-0.5).*A2.*(2+2*A+10*B));
%dedP_T2=dedr_T.*(rho2./Kt2);
dedT_P2=dedT_rho-dedr_T.*alpha2.*rho2;

A_phi2=1/3*(2*pi*Na*rho2).^.5.*(e.^2*(4*pi*eo*k*epsilon2.*tm2).^-1).^(3/2);
Ah2=-6*R*tm2.*A_phi2.*(1+tm2.*dedT_P2./epsilon2 +tm2.*alpha2/3);

% per paper, the D-H parameter for Apparent Cp is determined numerically
% since higher order derivatives of quantities are needed and difficult to
% code -  the finite difference seems adequately accurate.
Ac=(Ah2-Ah)/dT;

results.epsilon=reshape(epsilon,nP,nT);
results.dedP_T=reshape(dedP_T,nP,nT);
results.dedT_P=reshape(dedT_P,nP,nT);
results.Aphi=reshape(A_phi,nP,nT);
results.Ah=reshape(Ah,nP,nT);
results.Av=reshape(Av,nP,nT);
results.Ac=reshape(Ac,nP,nT);
results.Agam=reshape(A_phi*3,nP,nT);
results.rho=reshape(rho,nP,nT);
