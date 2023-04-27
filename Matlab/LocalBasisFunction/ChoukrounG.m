function results=ChoukrounG(PT,phase)
%Usage: [G,rho,Cp]=choukroun(PT,phase)
% PT is structure with P in MPa T in K
% phase is 'water' 'Ih' II III  V or VI
% output results has G (J/mol rho kg/m^3 and Cp J/K/kg

P=PT{1};
T=PT{2};

% Parameters from Choukroun and Grasset (2010)
% Below are the actual values in the program, x for values that differs from the paper
% This set of parameters do not reproduce the Vs reported in the paper for ices II and VI

%            Ih       II x    III x    liq       V        VI x
  data=[ 1.08600  0.85510  0.85460  0.81500  0.78300  0.74300
          273.160  248.850  256.430  400.000  273.310  356.150
         0.01900  0.02320  0.03750  0.10000  0.00500  0.02400
         0.00750  0.10800  0.02030  0.00500  0.01000  0.02000
         0.97400  0.99100  0.95100  1.00000  0.97700  0.96900
         0.03020  0.07930  0.09700  0.28400  0.12000  0.05000
         0.00395  0.00500  0.00102  0.00136  0.00160  0.00020] ;
% %            Ih       II x    III x    liq       V        VI x
%   V0   = [1.08600  0.85510  0.85460  0.81500  0.78300  0.74300] ;
%   Tref = [273.160  248.850  256.430  400.000  273.310  356.150] ;
%   a0   = [0.01900  0.02320  0.03750  0.10000  0.00500  0.02400] ;
%   a1   = [0.00750  0.10800  0.02030  0.00500  0.01000  0.02000] ;
%   b0   = [0.97400  0.99100  0.95100  1.00000  0.97700  0.96900] ;
%   b1   = [0.03020  0.07930  0.09700  0.28400  0.12000  0.05000] ;
%   b2   = [0.00395  0.00500  0.00102  0.00136  0.00160  0.00020] ;
% Below are the original values from the paper
data= [1.08600  0.84250  0.85500  0.81500  0.78300  0.74300
273.160  238.450  256.430  400.000  273.310  356.150
0.01900  0.06000  0.03750  0.10000  0.00500  0.02400
0.00750  0.00700  0.02030  0.00500  0.01000  0.00200
0.97400  0.97600  0.95100  1.00000  0.97700  0.96900
0.03020  0.04250  0.09700  0.28400  0.12000  0.05000
0.00395  0.00220  0.00200  0.00136  0.00160  0.00102];

switch phase
    case 'water'
      p=data(:,4); 
      p2=[0 0 0];
      c=[4190 9 -0.11];
      Tr=281.6;
    case 'Ih'
      p=data(:,1); 
      p2=[   251.16  209.9   -18.79];
      c=[74.11 7.56];
    case 'II'
      p=data(:,2); 
      p2=[252.32  300      -21.82];
      c=[2200 0];
    case 'III'
        p=data(:,3); 
        p2=[256.16  350.1    -16.19];
        c=[820 7];
    case 'V'
        p=data(:,5); 
        p2=[256.16  350.1    -17.43];
        c=[700 7.56];
    case 'VI'
        p=data(:,6); 
        p2=[273.31 632.4   -18.78];
        c=[940 5.5];
end
 
 To=p2(1);
 Po=p2(2);
 So=55.509*p2(3); %convert So to J/kg/K

 nT=length(T);
 nP=length(P);
 G=zeros(nP,nT);
 V=G;
 Cp=G;
 
 for i=1:nP
  for j=1:nT
     fV=@(x) 1e-3*(p(1)*(1+p(3)*tanh(p(4)*(T(j)-p(2))))*(p(5)+p(6)*(1-tanh(p(7)*x))));
     V(i,j)=fV(P(i));
     GV(i,j)=integral(fV,1e6*p2(2),1e6*P(i));
  end
 end
% 
if (~strcmp(phase,'water'))
    Cp=repmat((c(1)+c(2)*T(:)'),nP,1);
    Q1=repmat((c(1)*(T-To)+.5*c(2)*(T.^2-To^2)),nP,1);
    Q2=repmat((c(1)*log(T/To) + c(2)*(T-To)),nP,1);
else
    fCp=@(x) (c(1)+c(2)*exp(c(3)*(x-Tr)));
    fCpovT=@(x) (c(1)+c(2)*exp(c(3)*(x-Tr)))./x;
    Cp=repmat(fCp(T),nP,1);
    Q1=zeros(1,nT);
    Q2=Q1;
    for i=1:nT
        Q1(i)=integral(fCp,240,T(i));
        Q2(i)=integral(fCpovT,240,T(i));
    end
    Q1=repmat(Q1,nP,1);
    Q2=repmat(Q2,nP,1);
end
    G=(repmat(To*So,nP,nT)-So*repmat(T(:)',nP,1)+Q1-repmat(T(:)',nP,1).*Q2 +GV)/55.509; % eq 1 with conversion to J/mol
rho=V.^-1;

results.G=G;
results.Cp=Cp;
results.rho=rho;

