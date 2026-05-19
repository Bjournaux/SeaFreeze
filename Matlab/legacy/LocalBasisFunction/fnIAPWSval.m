function  out=fnIAPWSval(input,sp)
% function to return rho,vel, G, Cp, alpha K and Kp based on IAPWS  
%  ALL MKS with P in MPa
%    Usage: out=[rho,vel,G,Cp,alpha,S,U,H,K,Kp,mask]=fnIAPWSval(input)
%          where input is either a npts by (2 or 3) matrix of scatter points in [P  T] or input is a cell of {P,T}  
%                  P in MPa and T in K. 
%           rho in kg/m^3, vel in m/s, G in J/Mg Cp in J/kg/K alpha in K^(-1)
%           if sp is empty, it is loaded from file IAPWS_sp_strct.mat (needs to be in the path)
%           mask returns a matrix or vector with ones where this IAPWS spline is valid
%           and NaNs where not valid. A warning is issued for out of bounds
%           results.
%  Preloading the spline and passing it in the function call saves about 0.2 s or more rather 
%     than executing a load within this function
%
%  
%  calculations on a grid of PT values is much faster than scatter point calculations
%   
% JMB 2015  

mask=mk_mask_4_IAPWS(input);

%Determine whether gridded or scattered points
if iscell(input)                   % gridded data use b-form of spline
%   nd=length(input); 
   if(nargin==1),
       load IAPWS_sp_strct
       sp=sp_IAPWS;
   end   
    P=input{1};T=input{2};    
    %nP=length(P);
    %Tm=ones(nP,1)*T(:)';
    [Pm,Tm]=ndgrid(P,T);
    G=fnval(sp,{P,T});
    d1T=fnval(fnder(sp,[ 0 1]),{P,T});
    d2T=fnval(fnder(sp,[ 0 2]),{P,T});
    d1P=fnval(fnder(sp,[ 1 0]),{P,T});
    dPT=fnval(fnder(sp,[ 1 1]),{P,T});
    d2P=fnval(fnder(sp,[ 2 0]),{P,T});
    d3P=fnval(fnder(sp,[ 3 0]),{P,T});
else                               % scattered data in columns use pp-form of spline
   % [~,nd]=size(input);
    if(nargin==1),
       load IAPWS_sp_strct
       sp=fn2fm(sp_IAPWS,'pp');
    end  
    Tm=input(:,2);
    Pm=input(:,1);
    G=fnval(sp,input')';
    d1T=fnval(fnder(sp,[0 1]),input')';
    d2T=fnval(fnder(sp,[0 2]),input')';
    d1P=fnval(fnder(sp,[1 0]),input')';
    dPT=fnval(fnder(sp,[1 1]),input')';
    d2P=fnval(fnder(sp,[2 0]),input')';
    d3P=fnval(fnder(sp,[3 0]),input')';
end

out.S=-d1T;
out.rho=1e6*d1P.^(-1);  % 1e6 for MPa to Pa
out.G=G;
out.U=G-1e6*Pm./out.rho+Tm.*out.S;
out.H=out.U-Tm.*out.S;
out.Cp=-d2T.*Tm;
out.vel=real(sqrt(d1P.^2./(dPT.^2./d2T - d2P))); % MPa-Pa units conversion cancels
out.alpha=1e-6*dPT.*out.rho; % 1e6 for MPa to Pa
out.K=-d1P./d2P;
out.Kp=d1P.*d2P.^(-2).*d3P -1;
out.mask=mask;


function mask=mk_mask_4_IAPWS(PT)
if(iscell(PT))
    cell_flg=1;
    P=PT{1};
    T=PT{2};
else
    cell_flg=0;
    P=PT(:,1);
    T=PT(:,2);
end
PTmask=[300000 14000
    269000 14000
    237000 14000
    232000 13000
    222000 12400
    213000 11200
    204000 8700
    196000  7200
    180000  5400
    145000 4500
    97000  4300
    85000 4100
    73000 3000
    66000 2100
    62000 1600
    59000  1400
    32000   1000
    17000   820
    13000  650
    9800   600
    5400 500
    2500 375
    1400 303
    1000 270
    300 242
    273 242  %235
    0 242];

if cell_flg
    nP=length(P);
    nT=length(T);
    mask=ones(nP,nT);
    [Pm,Tm]=ndgrid(P,T);
    for i=1:nP
        Tc=interp1(PTmask(:,1),PTmask(:,2),P(i));
        id= T<Tc;
        mask(i,id)=nan;
    end
    id= Pm<50 & Tm>500;
    mask(id)=nan;
    id= Pm<100 & Tm>500 & Tm<900;
    mask(id)=nan;
else
    nT=length(T);
    mask=ones(nT,1);
    for i=1:nT
        Tc=interp1(PTmask(:,1),PTmask(:,2),P(i));
        if(T(i)<Tc),mask(i)=nan;end
        if(P(i)<50 && T(i)>500),mask(i)=nan;end
        if(P(i)<100 && T(i)>500 && T(i)<900),mask(i)=nan;end
    end
end
   
 id=find(isnan(mask), 1);
 if(not(isempty(id))),warning('some PT values are outside the valid range'),end
 
 