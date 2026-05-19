function Results=eqst(Data,Options)
% Usage:
%     Results=eqst(Data,Options)
%      
%       "Results" is a structure containing output values for rho, Cp, and G 
%       "Data" is a structure containing PT,rhoin,Cpin,vel, and Go
%          where: 
%               PT is a cell defining the P and T grid {P,T}
%               rhoin is a vector of densities at the reference isotherm, 
%               Cpin is a vector of specific heats at the referencef isotherm or is a grid of values over
%                 all pressures and temepratures - to force the Cp values
%               vel is a grid of sound speeds
%               Go is a  vector of Gibbs energy at the reference isotherm
%               (if not included Go is set to zero)
%
%       "Options" is a structure that contains all quantities associated with the
%           numerical analysis. If "Options" is not included in the input,
%           default values are used. "Options" can contain one or more elements
%           that changes default values. The list of values that can be set is
%          given below.
%             Tc:   defines the control points in T for spline fit to volumes vs temperature. 
%               default is 10 points spread over range of input T
%             nTc:  defines number of control points (default is 10)
%             logflg:  logflg=1 for logrithmic distribition of Tc (default is linear)
%             intflg:  default is trapezoidal integration of 1/velocities^2.  If
%               included in options, an interpolating cubic spline is integrated.
%             mask:  is a grid of ones or Nan's of size (nP,nT) to define the physical regime
%            loglam_start:  log10(lam) for the initial damping parameter (default=1)     
%            loglam_end:   log10(lam) for the damping parameter at the highest pressures (default=3)    
%            pflg:    if pflg=1 plots quantities associated with the adiabatic to isothermal
%                correction - provides information if the integration is not working
%            mdrv:  derivative of V vs T that is to be smooth (default=4)
%            dirctn_flg:  =1 for integrating up in pressure and not 1 for integrating down in pressure
%
%      This algorithm can suffer from the growth of numerical instabilities
%      at higher pressures.  Two stabilizing ideas are incorporated. 
%            (1) With each step in pressure, the damping of V(T) fits is
%            increased by a logrithmic increment defined by the range between loglam_start and loglam_end.
%            (2) The second derivatives of V(T) are averaged over 5 points 
%                to reduce growth of point to point instabilities (may be
%                commented out)
%       Setting  "pflg" to 1 provides plots at each pressure step.  This
%       highlights at what pressure and over what temperature range the
%       algorithm fails. This can then provide guidance in modification of
%       the input sound speed surface.  A defect in the sound speed surface 
%       is a typical source of integration difficulty.
%       
% use MKS units (P in MPa)
%  JMB 2014 - 2017

if nargin==5
    Options=[];
end

if isfield(Options,'mask')
    mask=Options.mask;
else
    mask=ones(size(Data.vel));
end

if isfield(Options,'loglam_start')
    loglam_start=1+Options.loglam_start;
else
    loglam_start=1;
end

if isfield(Options,'loglam_end')
    loglam_end=1+Options.loglam_end;
else
    loglam_end=3;
end

if isfield(Options,'pflg')
    pflg=1;
else
    pflg=0;
end

if isfield(Options,'mdrv')
    mdrv=Options.mdrv;
else
    mdrv=4;
end

if isfield(Options,'ordr')
    ordr=Options.ordr;
else
    ordr=6;
end
if isfield(Options,'nReg')
    nReg=Options.nReg;
else
    nReg=1;
end

if isfield(Options,'nTc')
    nTc=Options.nTc;
else
    nTc=10;
end

if isfield(Options,'logflg')
    logflg=1;
else
    logflg=0;
end

if isfield(Options,'Tc')
    Tc=Options.Tc;
    nTc=length(Tc);
end

if isfield(Options,'intflg')
    intflg=Options.intflg;
else
   intflg='numerical';
end

if isfield(Options,'Cpsmoothflg')
    Cpsmoothflg=1;
else
    Cpsmoothflg=0;
end

if isfield(Options,'dirctn_flg')
    dirctn_flg=Options.dirctn_flg;
else
   dirctn_flg=1;
end

kT=ordr;
P=Data.PT{1};
T=Data.PT{2};
T=T(:)'; % enforce T as row
nT=length(T);
facP=1e6;  %conversion of P from MPa to Pa;
P=facP*P(:); % insist that P is a column and in Pa 
nP=length(P);
dloglam=(loglam_end-loglam_start)/nP;
lam=10^loglam_start;

if dirctn_flg==1
    id_start=1;
else
    id_start=nP;
end

if logflg
    Tc=logspace(log10(min(T)-1),log10(max(T)+1),nTc);
else
    Tc=linspace(min(T)-1,max(T)+1,nTc);
end

if(nTc>.5*nT),error('Too many control points relative to temperature data points'),end

velsq= mask.*(Data.vel).^(-2);
velsq(isnan(mask))=0;
switch intflg
    case 'splinefit'
        spv=csapi({P, T},velsq); % assumption that vel is defined everywhere on the grid in P and T. 
        %                           csapi is a MATLAB provided cubic spline interpolations function
        velsqInt=fnval(fnder(spv,[-1 0]),{P,T}); % analytically evaluate the integral of the spline representing 1/vel^2
    case 'numerical'
        velsqInt=cumtrapz(P,velsq);  % this numerical integration should handle the situation of vel not defined everywhere
end

dvelsqInt=diff(velsqInt);
Tr=mkgrid(Tc,nReg);
nTr=length(Tr);
nv=min(size(Data.Cpin));
if(nv==1)
    fit_Cp=1;
else
    fit_Cp=0;
end

wtR=area_wt(Tr);

% set up b spline for V(T)
Tknts=optknt(Tc,ordr);
Tcol = spcol(Tknts,ordr,brk2knt(T,ordr));
Vcol=Tcol(1:ordr:end,:);
%d2Vcol=Tcol(3:ordr:end,:);
Tcolreg = spcol(Tknts,ordr,brk2knt(Tr,ordr));
RegCol=Tcolreg((mdrv+1):ordr:end,:);  % specified derivative of V with T to be made small
NA=norm(Vcol,1);
NR=norm(RegCol,1);

sp.form='B-';
sp.knots=Tknts;
sp.number=length(Tc);
sp.order=kT;
sp.dim=1;


rho=zeros(nP,nT);
alpha=rho;
Cp=rho;
G=rho;
cor=rho;
dVdT=rho;
d2VdT2=rho;

%set fit flag for specific heat 
if fit_Cp
        Cp(id_start,:)=Data.Cpin(:)';
else
    Cp=Data.Cpin; 
end

rho(id_start,:)=Data.rhoin(:)';
if isfield(Data,'Go')
   G(id_start,:)=Data.Go(:)';
end
    
V=rho(id_start,:).^(-1);
Vnorm=std(V)^-1;

A=[Vnorm*Vcol;lam*NA/NR*(wtR(:)*ones(1,nTc)).*RegCol];
B=[Vnorm*V(:);zeros(nTr,1)];
sp.coefs=(A\B)'; % solve for V(T)

dVdT(id_start,:)=fnval(fnder(sp),T);
d2VdT2(id_start,:)=fnval(fnder(sp,2),T);
alpha(id_start,:)=dVdT(id_start,:).*rho(id_start,:);
cor(id_start,:)=alpha(id_start,:).^2./Cp(id_start,:).*T;
incld=find(not(isnan(mask(id_start,:))));

% plot quantities if debugging
if pflg
    figure(100)
    txt=sprintf('P= %6.0f ',P(id_start)/1e6);
    subplot(411)
    plot(T(incld),rho(id_start,incld))
      title('Density')
    rm=max(rho(id_start,incld));
    text(300,rm-10,txt)
     subplot(412)
    plot(T(incld),Cp(id_start,incld),'-o')
      title('Cp')
    subplot(413)
    plot(T(incld),d2VdT2(id_start,incld),'-o')
      title('d2VdT2')
    subplot(414)
    plot(T(incld),cor(id_start,incld))
    title('Correction')
    pause
end

% begin integration to high pressure - starting at the second pressure
% point

if dirctn_flg==1
    i_index=2:nP;
    cor_fac=.8;
else
    i_index=(nP-1):-1:1;
    cor_fac=1.2;
end

for i=i_index
    icount=0;
    test=1;
    incld=find(not(isnan(mask(i,:))));
%    lam=lam_end*lam;
    lam=lam*10^dloglam;
% increment densities using velocity integral and adiabatic correction term
    if dirctn_flg==1
        cor(i,incld)=cor_fac*cor(i-1,incld);   % predict that the correction for next step is smaller than previous step 
    else
        cor(i,incld)=cor_fac*cor(i+1,incld);
    end
while (test>1e-14 && icount<20) % stop when corrections do not change from iteration to iteration
    icount=icount+1;
    if dirctn_flg==1
       rho(i,incld)= rho(i-1,incld)+dvelsqInt(i-1,incld)+trapz(P(i-1:i),cor(i-1:i,incld)); % haven't found a stable way to remove this trapezoidal integration  
    else
        rho(i,incld)= rho(i+1,incld)-dvelsqInt(i,incld)+trapz(P(i+1:-1:i),cor(i+1:-1:i,incld)); % haven't found a stable way to remove this trapezoidal integration
    end
    V=rho(i,incld).^(-1);
    Vnorm=std(V)^-1; 
    A=[Vnorm*Vcol(incld,:); lam*NA/NR*(wtR(:)*ones(1,nTc)).*RegCol];
    B=[Vnorm*V(:); zeros(nTr,1)];   
    sp.coefs=(A\B)';
    dVdT(i,incld)=fnval(fnder(sp),T(incld));
    d2VdT2(i,incld)=fnval(fnder(sp,2),T(incld));   
    alpha(i,(incld))=dVdT(i,(incld)).*rho(i,(incld));
    
    if fit_Cp  % section to updata Cp
    %following does trapezoidal integration for Cp at new pressure
        if dirctn_flg==1
            Cp(i,incld)=Cp(i-1,incld) - trapz(P(i-1:i),(ones(2,1)*T(incld)).*d2VdT2(i-1:i,incld));
        else
            Cp(i,incld)=Cp(i+1,incld) - trapz(P(i+1:-1:i),(ones(2,1)*T(incld)).*d2VdT2(i+1:-1:i,incld)); 
        end  
    
  % following is a 5 point averaging of Cp vs T to damp high frequency numerical instabilities
  if Cpsmoothflg
         for j=3:(length(incld)-2)
             Cp(i,incld(j))=mean(Cp(i,incld(j-2:j+2)));
         end  
  end

% Calculate new value of correction term
       corN=alpha(i,incld).^2./Cp(i,incld).*T(incld);
       test=max(abs(cor(i,incld)-corN));
       cor(i,incld)=corN;
    
% plot stuff if debugging:    
  if pflg
    figure(100)
    txt=sprintf('P= %6.0f ',P(i)/1e6);
    subplot(411)
    plot(T(incld),rho(i,incld))
      title('Density')
    rm=max(rho(i,incld));
    text(300,rm-10,txt)
     subplot(412)
    plot(T(incld),Cp(i,incld),'-o')
      title('Cp')
    subplot(413)
    plot(T(incld),d2VdT2(i,incld),'-o')
      title('d2VdT2')
    subplot(414)
    plot(T(incld),trapz(P(i-1:i),cor(i-1:i,incld)),T(incld), dvelsqInt(i-1,incld),'--')
    title('Correction')
    pause
  end   % end of debug plotting
  
    end   %end of updating Cp

    if dirctn_flg==1
        G(i,incld)=G(i-1,incld)+ trapz(P(i-1:i),rho(i-1:i,incld).^-1);
    else
        G(i,incld)=G(i+1,incld)+ trapz(P(i+1:-1:i),rho(i+1:-1:i,incld).^-1);
    end
end  % end of iterations on correction term
end % end of pressures loop

Results.rho=rho;
Results.Cp=Cp;
Results.G=G;

 
 function wt=area_wt(PT)
% this subfunction assigns a weight to every regularization point based on
% the volume encompassed by that point
if(iscell(PT))
 if (length(PT)==2)   
    p=PT{1};
    dp=diff(p(:))';
    t=PT{2};
    dt=diff(t(:))';
    wtp=.5*[dp(1) dp(1:end-1)+dp(2:end) dp(end)];
    wtt=.5*[dt(1) dt(1:end-1)+dt(2:end) dt(end)];
    wt=wtp(:)*wtt;
    wtsum=sum(wt(:));
    wt=sqrt(wt/wtsum);
 elseif( length(PT)==3)  
      p=PT{1};
      dp=diff(p(:))';
      t=PT{2};
      dt=diff(t(:))';
      m=PT{3};
      dm=diff(m(:))';
    wtp=.5*[dp(1) dp(1:end-1)+dp(2:end) dp(end)];
    wtt=.5*[dt(1) dt(1:end-1)+dt(2:end) dt(end)];
    wtm=.5*[dm(1) dm(1:end-1)+dm(2:end) dm(end)];
    [a,b,c]=ndgrid(wtp,wtt,wtm);
    wt=a.*b.*c;
    wtsum=sum(wt(:));
    wt=sqrt(wt/wtsum);
  end
 else    
    dx=diff(PT(:))';
    wtx=.5*[dx(1) dx(1:end-1)+dx(2:end) dx(end)];
    wtsum=sum(wtx(:));
    wt=sqrt(wtx/wtsum);
end
 
function xg=mkgrid(x,nreg)
% this subfunction puts nreg points between every value of x
xg=x(:)';
dx=diff(xg);
if nreg>1
for i=2:nreg
    xg=[xg x(1:end-1)+(i-1)*dx/nreg];
end
end
xg=sort(xg);

