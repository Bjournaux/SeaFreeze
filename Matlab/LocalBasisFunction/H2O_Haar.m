function results= H2O_Haar(input)
% use Haar function to get useful thermodynamic properties
% Haar_out=H2O_Haar(PT)
% PT can be structure of {P,T} or matrix of [P(:) T(:)]
% Input in MPa and K
% output in standard units
% note that the mex file has G based on Berman 1988.  This code re-centers
% energies and entropies such that at the triple point (T=273.16 P = 611.655 Pa)
% U (internal energy) and S (specific entropy) are zero.
% JMB 2017

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

% some units and energy standards
    P_bar=10*pm(:)';
    T_K=tm(:)';
    mw = 0.0180152;  %kg/mol
    % triple point values to recenter EOS
    ho=.611872;  % IAPWS value for H at triple point
    %pr=611.655e-6;
    %Tr=273.16; 
    Sr=63.338204089922485;% value from HAAR.c
    Hr=-2.877188194348947e+05;  % value from HAAR.c

% call the mex file:
    [p,~,H,S,Cp,~,v,dvdt,dvdp]=Haar(P_bar,T_K);

% convert to our units
    rhoH=1e6*mw*v(:).^(-1);
    P2=p*1e5;
    Cp=Cp(:);
    alpha=dvdt(:)./v(:);
    Kt=-dvdp(:).^(-1).*v(:)*1e5/1e9;  % corrected for bars and for GPa
    Cv=Cp(:)-T_K(:).*(v(:)*1e-6).*alpha(:).^2.*(Kt(:)*1e9);
    Ks=Kt(:).*Cp(:)./Cv(:);
    gamma=alpha.*Kt.*(v(:)*1e3)./Cv;
    velocity=sqrt(1e9*Ks(:)./rhoH(:));
    Cp=Cp/mw;
    Cv=Cv/mw;

%Re-center energies
    H=(H/mw-Hr/mw+ho);   
    S=(S-Sr)/mw;
    G=H-T_K'.*S;

% put results in structure
    results.rho=reshape(rhoH,nP,nT);
    results.alpha=reshape(alpha,nP,nT);
    results.G=reshape(G,nP,nT);
    results.S=reshape(S,nP,nT);
    results.vel=reshape(velocity,nP,nT);
    results.Cp=reshape(Cp,nP,nT);
    results.Cv=reshape(Cv,nP,nT);
    results.gamma=reshape(gamma,nP,nT);
    results.Ks=reshape(Ks,nP,nT);
    results.Kt=reshape(Kt,nP,nT);
    results.H=reshape(H,nP,nT);
    results.P=reshape(P2/1e6,nP,nT);












