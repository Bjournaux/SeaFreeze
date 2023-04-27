function pabalan=PabalanEOS(input,pflg,H2Oflg)
%this returns the Na2SO4 EOS based on Pabalan and Pitzer
% Geochimica et Cosmochimica ACla Vol. S2. pp. 2393-2404
% Heat capacity and other thermodynamic properties of Na2SO4(aq) in hydrothermal solutions
% and the solubilities of sodium sulfate minerals in the system Na-CI-SO .. -OH-H20 to 300°C
%  Usage:  pabalan=PabalanEOS(input,pflg,H2Oflg)
%        where the output structure contains standard EOS quantities and   activities
%        including Gs (J/kg) Gibbs energy per kg of solution and G per kg of water
%       input:   {T,m} or [T(:) m(:)] with T in (K) and concentration in kg/mol pflg=1 for 200 bars and 0 for 1 bar
%        pflg:  1 for 20 MPa and 0 for .1 MPa - the two pressure reported in the study.
%        H2Oflg 'Haar' or 'IAPWS'  the first uses Bradley and Pitzer for
%        the dielectric constant the second uses Fernandez et al 1997
%  Values still need to be checked against an absolute G based on solubility
% JMB 2017

% organize input
    if(iscell(input))
        T=input{1};
        m=eps+input{2};
        [tm,mm]=ndgrid(T,m);
        nT=length(T);
        nm=length(m);
    else
        tm=input(:,1);
        tm=input(:,1);
        mm=eps+input(:,2);
        nT=length(tm);
        nm=1;
    end
 tm=tm(:);
 mm=mm(:);
 npts=length(tm);

% some constants:
b=1.2;
a1=1.4;
R=8.3144621;
MW=0.14204;
omega=1000/18.0152;
Tr=298;
fac=(1+MW*mm);  % conversion factor for per kilo of water and per kilo of solution

I=3*mm;
sqrtI=sqrt(I);
x=a1*sqrtI;
g1=2*x.^-2.*(1-(1+x).*exp(-x));

if pflg
    P=20;
else
    P=.1;
end
switch H2Oflg
    case 'Haar'
        results = Bradley1979dielectric({P,T});
        Cpw=results.haar.Cp(:)*ones(1,nm);
        muo_w=results.haar.G(:)*ones(1,nm); % in J per kg of water
    case 'IAPWS'
        results = Fernandez1997dielectric({P,T});
        muo_w=results.iapws.G(:)*ones(1,nm);
        Cpw=results.iapws.Cp(:)*ones(1,nm);
end        
muo_w=muo_w(:);
Cpw=Cpw(:);
Aphi=results.Aphi(:)*ones(1,nm);
Aphi=Aphi(:);
Av=results.Av(:)*ones(1,nm);
Av=Av(:);
Ah=results.Ah(:)*ones(1,nm);
Ah=Ah(:);
Ac=results.Ac(:)*ones(1,nm);
Ac=Ac(:);
tm=tm'; % this transpose is temporary and makes the basis functions easier to read
basis_G=[  (tm.^2)/6 
            tm/2 
            tm.^2.*(log(tm)-5/6)/6 
            (tm.^3)/12 
            (tm.^4)/20 
            (tm/2+3*(227^2)/2*(tm.^-1)+227*(tm-227).*log(tm-227).*(tm.^-1)) 
            -(tm/2 +3*(647^2)/2.*(tm.^-1)-647*(647-tm).*log(647-tm).*(tm.^-1))
            -tm.^-1        
            -(Tr^2)*(tm.^-1)
             ones(1,npts) 
             ones(1,npts)
]';
         
 basis_H=[  tm/3 
            ones(1,npts)/2 
            tm.*(log(tm)-1/3)/3 
            (tm.^2)/4 
            (tm.^3)/5 
            tm.^-2.*(((tm-227).^2)/2 +454*(tm-227) +227^2*log(tm-227))      
            tm.^-2.*(-(647-tm).^2/2+1294*(647-tm)-647^2*log(647-tm))
            tm.^-2
            (Tr^2)*(tm.^-2)
            zeros(1,npts)
            zeros(1,npts)
            ]';

basis_Cp=[  ones(1,npts)
            tm.^-1
            log(tm)
            tm
            tm.^2
            (tm-227).^-1
            (647-tm).^-1
            zeros(1,npts)
            zeros(1,npts)
            zeros(1,npts)
            zeros(1,npts)
      ]';           
 tm=tm';   %putting tm back to original orientation   

parm_1bar=[
2.1549644e-2      3.6439508          -4.6590760e-2
-7.6918219e-1    -9.1962646e1         1.1711403
-3.5486084e-3    -6.5961772e-1        8.4319701e-3
4.0811837e-6      1.6043868e-3        -2.0439550e-5
0.0              -6.5972836e-7         8.3348147e-9
0.0               1.6491982e-1         -1.3141008e-3
0.0               2.2057312e-1         -1.6791063e-3
3915.434531       708447.986899        -6717.929066
1.738512e-3       5.820066e-3          -1.117462e-4
55.769488          5768.102375         -68.202263
-1.255087e-2      7.037660e-1          3.808550e-3
];


data=[
 1.7278418e6     6.0955633e-1    1.1040235     -1.0300087e-1
 -4.6433635e7    -1.6090797e1    -2.5758534e1    2.6412809
 -3.0786428e5    -1.0932828e-1    -2.0290775e-1  1.8580153e-2
 6.8028314e2      2.5321479e-4     5.3309441e-4   -4.4349928e-5
 -2.4386675e-1    -9.9384034e-8     -2.3576724e-7  1.7883333e-8
 1.0033867e5      4.0107638e-2      0.0           -5.1836594e-3
 -2.5274859e5      2.1711348e-2     1.4455381e-1   -3.8805720e-3
     0            92308.895357    363078.716679   -15447.360134
  0                1.722469e-3      5.512612e-3    -1.152988e-4
      0             963.974106     1926.602872      -155.941254
        0                 -1.04936e-2     6.90077e-1      3.74436e-3
];

parm_200bar=data(:,2:4);
if pflg
    mod=basis_G*parm_200bar;
    mod_Cp=basis_Cp*parm_200bar;
else
    mod=basis_G*parm_1bar;
    mod_Cp=basis_Cp*parm_1bar;
end
Cpo=basis_Cp*data(:,1);
Bmxo=mod(:,1);
Bmx1=mod(:,2);
Cmx=mod(:,3);
Bmxoj=mod_Cp(:,1);
Bmx1j=mod_Cp(:,2);
Cmxj=mod_Cp(:,3);

Bmx_phi=Bmxo+Bmx1.*exp(-x);
Bmx= Bmxo  +   Bmx1.*g1;
Bmx_g=Bmx+Bmx_phi;
Bmx_j=Bmxoj+Bmx1j.*g1;
Cmx_phi=Cmx;

% here is the calculation of the excess and apparent quantities
phi=-2*Aphi.*sqrtI./(1+b*sqrtI)+4/3*mm.*Bmx_phi +16/3*mm.^2.*Cmx +1;
gam=exp(-2*Aphi.*(sqrtI./(1+b*sqrtI)+2/b*log(1+b*sqrtI))+4/3*mm.*Bmx_g +8*mm.^2.*Cmx);
phi_Cp= Cpo + 3*Ac.*log(1+b*sqrtI)/b - 4*R*tm.^2.*(mm.*Bmx_j+2*mm.^2.*Cmxj);
% Cp1= 3*Ac.*log(1+b*sqrtI)/b ;
% Cp2=- 4*R*tm.^2.*(mm.*Bmx_j+2*mm.^2.*Cmxj);

Gex=R*(tm).*((-4*(Aphi.*I).*(log(1+b*sqrtI)))/b + 2*2*((mm.^2).*Bmx + 2*(mm.^3).*(Cmx)));
% G1=R*(tm).*((-4*(Aphi.*I).*(log(1+b*sqrtI)))/b);
% G2=R*(tm).*( 2*2*((mm.^2).*Bmx + 2*(mm.^3).*(Cmx)));
%Gex2=3*R*tm.*mm.*(log(gam)-(phi-1)); % redundancy check - it agrees to 10-12 so I think the Gex terms are correct
% max(max(Gex-Gex2))
% calculate solution Cp (J/K/ kg of solution)

Cpw=Cpw(:);
Cps=(Cpw+phi_Cp.*mm)./fac;

% here is the "non-excess" compositional dependent part of G
Gnon=R*3*(tm.*mm).*((log(mm)) -1);  %+log(2)/3

% need standard state properties of salt and water

% Here is calculatiom for salt based on Pabalan's integration constants
muoos=R*298.15*(-511.414); % energy per mole at 1 bar and 298.15 K
Sos=R*16.263;  % entropy at 1 bar and 298.15 K
% next integrate standard state Cp and Cp/T to get energy at other
% temperatures
Q1=zeros(npts,1);
Q2=Q1;
f2 = @(x) (PabalanCp(x)'./x) ; % ACp/T
f1= @(x) PabalanCp(x)';  % ACp
for i=1:npts  
        Q1(i)=integral(f1,298.15,tm(i));  % integrate from 298 to T
        Q2(i) = integral(f2,298.15,tm(i));
end

muo_s=1.27e6+(muoos-Sos*(tm-298.15)+Q1-tm.*Q2); % this takes a value at 298 and 1 bar and 
mu_s=muo_s + R*3*tm.*log(gam.*mm);  % arbitrary offset of 1.27e6 to slide chemical potential to near zero

lnaw=-mm.*phi*3/omega;
mu_w=muo_w/omega+R*tm.*lnaw;
% this alternate definition has been tested - it works
%G2=omega*mu_w + mm.*mu_s;  % this needs to be equal to the G defined below

% now add them up
Go= muo_w + mm.*muo_s;
G=Gex+Go+Gnon;  % per kg of water
Gs=G./fac;   % per kg of solution

% put in otput structure
pabalan.G=reshape(G,nT,nm);
%pabalan.G2=reshape(G2,nT,nm);
pabalan.Gs=reshape(Gs,nT,nm);
pabalan.Gex=reshape(Gex,nT,nm);
pabalan.mu_w=reshape(mu_w,nT,nm);
pabalan.muo_w=reshape(muo_w,nT,nm);
pabalan.mu_s=reshape(mu_s,nT,nm);
pabalan.muo_s=reshape(muo_s,nT,nm);
pabalan.Gnon=reshape(Gnon,nT,nm);
pabalan.phi=reshape(phi,nT,nm);
pabalan.gamma=reshape(gam,nT,nm);
pabalan.phi_Cp=reshape(phi_Cp,nT,nm);
pabalan.Cpo=reshape(Cpo,nT,nm);
pabalan.Cps=reshape(Cps,nT,nm);
pabalan.fac=reshape(fac,nT,nm);
% pabalan.Cp1=reshape(Cp1,nT,nm);
% pabalan.Cp2=reshape(Cp2,nT,nm);
pabalan.Cpo=reshape(Cpo,nT,nm);
% pabalan.G1=reshape(G1,nT,nm);
% pabalan.G2=reshape(G2,nT,nm);

function Cpo=PabalanCp(T)
% input T in K
T=T(:);
nt=length(T);
basis=[ones(nt,1) T.^-1 log(T) T T.^2 (T-227).^-1 (647-T).^-1];
pCpo=[
 1.7278418e6 
-4.6433635e7 
-3.0786428e5 
6.8028314e2
-2.4386675e-1
 1.0033867e+5 
-2.5274859e5 
];
Cpo=basis*pCpo;

