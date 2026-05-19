function spG=spgmft(PTm,G,rho,Cp,options)
% create a spline for Gibbs energy - use G (optionally rho and Cp) as constraint and
% regularize in the non-physical region with smooth mdrv th derivatives. 
%    G can be function of P and T or function of P T and composition (m)
%
%    Units for all quantities are MKS!! (J/kg, ,MPa, K kg/m^3)  
%      (internal factor of 10^6  converts  MPa input to Pa)
% 
%    G rho and Cp are grids of values with PTm containing cells with P and T values.  
%        Alternative is to enter G along a base -
%         In that case G is a vector of values and Cp and rho are used to reconstruct G by colocation at higher pressures
%
%    The spline is on a rectangular grid of points. In the real world, 
%    phase boundaries interfere and sections of the grid represent metastable regions where 
%    no reasonable thermodynamic values are available.  Thus, to solve the problem on the grid (insure that the
%    matrixes are adequately conditions and not too rank-deficit) regularization is applied 
%        - assert that the surfaces need to be "smooth" where data are not available. 
%
% function call:
%      sp_G=spgft(PTm,G,rho,Cp,options)
%    input "PTm" is a cell containing P and T or P T and m
%    input "G" is either gridded data at PTm values or Tm values for P=Pr
%    inputs "rho" and "Cp" are either gridded data at PTm values or Tm values for P=Pr
%
%    all remaining inputs are in a structure called "options" and have default values if not provided in the structure
%       "mask" has ones where input thermodynamic properties are valid and NaN where they are not (Default is all ones)
%       "PTmc" is a cell containing the control points for the spline (Default is the PTm grid)              
%       "lam" is vector of smoothing parameters (default is all ones )
%       "weight" is vector of weights for the relative importance of densities and Cp relative to G (default is [1 1]);
%       "mdrv" is vector of the derivatives to use in smoothing (default is [4 4 3]
%       "ordr" is vector of integer order's for the spline  (default is [6 6 4])
%       "nReg" sets the number of regularization points for each control point (default is  1)
%       "normflg" =1 to rescale data based on their standard deviations (default is nop scaling)
%       "algorithm" to use Least Square fitting (default is Normal equations)
%
% custom functions mkCmn wt_area and mkgrid are nested in this file
%  JMB 2015

ndim=length(PTm);

m_flg=0;
if ndim==3
    m_flg=1;
end

if m_flg
    if length(size(rho))~=3
        error('density array is wrong size')
    end
    if length(size(Cp))~=3
        error('Cp array is wrong size')
    end
else
    if length(size(rho))~=2
        error('density array is wrong size')
    end
    if length(size(Cp))~=2
        error('Cp array is wrong size')
    end
end

% extract P and T from cell, convert P to Pa units
P=PTm{1}*1e6;
P=P(:);
T=PTm{2};
T=T(:)';
if m_flg
    m=PTm{3};
    m=m(:)';
end

if(isfield(options,'PTmc'))
    PTmc=options.PTmc;
    Pc=PTmc{1}*1e6;
    Tc=PTmc{2};
    if m_flg
        mc=PTmc{3};
    end
else
    Pc=P;
    Tc=T;
    if m_flg
        mc=m;
    end    
end

if(isfield(options,'normflg'))
    normflg=1;
else
    normflg=0;
end

reg_flg=1;
if(isfield(options,'nReg'))
    nReg=options.nReg;
    if(isempty(nReg))
      reg_flg=0;
    end
else
    nReg=ones(1,ndim);
end

if(isfield(options,'ordr'))
    ordr=options.ordr;
    kp=ordr(1);
    kt=ordr(2);
    if m_flg
        km=ordr(3);
    end
else
    if ndim==2
      %ordr=[6 6];
       kt=6;
       kp=6;
    elseif ndim==3
      %ordr=[6 6 4];
      kt=6;
      kp=6;
      km=4;
    end

end

if(isfield(options,'mdrv'))
    mdrv=options.mdrv;
else
    if ndim==2
      mdrv=4*ones(1,2);
    elseif ndim==3
        mdrv=[4 4 3];
    end
end

if(isfield(options,'lam'))
    lam=options.lam;
else
    if ndim==2
      lam=[1 1];
    elseif ndim==3
        lam=[1 1 1];
    end
end

if(isfield(options,'mask'))
    mask=options.mask;
else
    mask=ones(size(rho));
end

if(isfield(options,'weight'))
    weight=options.weight;
else
    weight=[1 1];
end

if(isfield(options,'algorithm'))
    algorthm='LeastSqr';
else
    algorthm='NormlEqn';
end

% determine number of data and control points
nP=length(P);
%nT=length(T);
nPc=length(Pc);
nTc=length(Tc);
if m_flg
    nm=length(m);
    nmc=length(mc);
end


mdrv=mdrv+1; %increment derivatives to align chosen regularization derivative with the order of basis functions returned by spcol

% find indexes for regions where G is defined and regions where it is not
id_incl=find(not(isnan(mask(:))));

% determine whether density or Cp data are included
if(isempty(rho)) 
    rhoflg=0;
else
    rhoflg=1;
end

if(isempty(Cp))
    Cpflg=0;
else
    Cpflg=1;
end

% determine if G is defined only at the base or everywhere
if m_flg
   if (length(size(squeeze(G))))==2
        Gvflg=1;
   else
        Gvflg=0;
   end 
else
    if min(size(squeeze(G)))==1
        Gvflg=1;
    else
        Gvflg=0;
    end
end

m_rb=mask(1,:,:);
id_rb=find(not(isnan(m_rb(:))));
n_dat=length(id_incl);
norm_fac=sqrt(n_dat);  % will weight data vs regularization by square root of number of data points

if Gvflg % handle G only at 1 bar or G for all P and T and m
    G_norm=std(G(id_rb))^-1; % use standard deviation of G at 1 bar for normalization
else
    G_norm=std(G(id_incl(:)))^-1; % use standard deviation of G(P,T) for normalization
end

% rho and Cp expressed as linear derivatives of G and determine standard
% deviations of the "data" derivatives
if rhoflg, Gp=rho.^-1;Gp_norm=std(Gp(id_incl(:)))^-1;end    % use standard deviation of Gp for normalization
if Cpflg
    if m_flg
        Gtt=-Cp./repmat(T(:)',nP,1,nm);Gtt_norm=std(Gtt(id_incl(:)))^-1;
    else
        Gtt=-Cp./repmat(T(:)',nP,1);Gtt_norm=std(Gtt(id_incl(:)))^-1;   
    end   % use standard deviation of Gtt for normalization
end

Pknts=aptknt(Pc(:)',kp);
Tknts=aptknt(Tc(:)',kt);
Tcol = spcol(Tknts,kt,brk2knt(T,kt));
Pcol = spcol(Pknts,kp,brk2knt(P,kp));
if m_flg
    mknts=aptknt(mc(:)',km);
    mcol = spcol(mknts,km,brk2knt(m,km));
end

% create the matrixes of basis functions
% depends on whether G is defined everywhere or just at the base
if Gvflg
    if m_flg
        RdatG=makeCmn(Pcol(1,:),Tcol(1:kt:end,:),mcol(1:km:end,:));
    else
        RdatG=makeCmn(Pcol(1,:),Tcol(1:kt:end,:));
    end
else
    if m_flg
        RdatG=makeCmn(Pcol(1:kp:end,:),Tcol(1:kt:end,:),mcol(1:km:end,:));
    else
        RdatG=makeCmn(Pcol(1:kp:end,:),Tcol(1:kt:end,:));
    end
end
% basis functions for densities (volumes) and specific heats:
if m_flg
     RdatV=makeCmn(Pcol(2:kp:end,:),Tcol(1:kt:end,:),mcol(1:km:end,:));
     RdatC=makeCmn(Pcol(1:kp:end,:),Tcol(3:kt:end,:),mcol(1:km:end,:));
else
    RdatV=makeCmn(Pcol(2:kp:end,:),Tcol(1:kt:end,:));
    RdatC=makeCmn(Pcol(1:kp:end,:),Tcol(3:kt:end,:));
end

% in the case that regularization is used, set up the basis functions
% then weight by volume contained around each regularization point
if reg_flg
    Pr=mkgrid(Pc,nReg(1));
    Tr=mkgrid(Tc,nReg(2));
    Tcolr = spcol(Tknts,kt,brk2knt(Tr,kt));
    Pcolr = spcol(Pknts,kp,brk2knt(Pr,kp));
    if m_flg
        mr=mkgrid(mc,nReg(3));
        mcolr = spcol(mknts,km,brk2knt(mr,km));
        wtreg=area_wt({Pr,Tr,mr});
        RregT=makeCmn(Pcolr(1:kp:end,:),Tcolr(mdrv(2):kt:end,:),mcolr(1:km:end,:));    
        RregP=makeCmn(Pcolr(mdrv(1):kp:end,:),Tcolr(1:kt:end,:),mcolr(1:km:end,:));   
        Rregm=makeCmn(Pcolr(1:kp:end,:),Tcolr(1:kt:end,:),mcolr(mdrv(3):km:end,:));  
        [nr,nc]=size(RregT);
        [it,jt,xt]=find(RregT);    
        [ip,jp,xp]=find(RregP);
        [im,jm,xm]=find(Rregm);
        RregT=sparse(it,jt,wtreg(it).*xt,nr,nc);      
        RregP=sparse(ip,jp,wtreg(ip).*xp,nr,nc);
        Rregm=sparse(im,jm,wtreg(im).*xm,nr,nc);
    else
        wtreg=area_wt({Pr,Tr});
        RregT=makeCmn(Pcolr(1:kp:end,:),Tcolr(mdrv(2):kt:end,:));    
        RregP=makeCmn(Pcolr(mdrv(1):kp:end,:),Tcolr(1:kt:end,:));   
        [nr,nc]=size(RregT);
        [it,jt,xt]=find(RregT);    
        [ip,jp,xp]=find(RregP);
        RregT=sparse(it,jt,wtreg(it).*xt,nr,nc);      
        RregP=sparse(ip,jp,wtreg(ip).*xp,nr,nc);
    end
end

% alternate way to normalize - does not seem to help
% NG=norm(RdatG,1);
% NV=norm(RdatV,1);
% NC=norm(RdatC,1);
% NRP=norm(RregP,1);
% NRT=norm(RregT,1);
% NRm=norm(Rregm,1);
% Gp_norm=Gp_norm*NG/NV;
% Gtt_norm=Gtt_norm*NG/NC;
% norm_facP=NG/NRP;
% norm_facT=NG/NRT;
% norm_facm=NG/NRm;

norm_facP=norm_fac;
norm_facT=norm_fac;
norm_facm=norm_fac;

% set up the spline structure for the output
spG.form='B-';
if m_flg
    spG.knots={Pknts,Tknts,mknts};
    spG.number=[nPc nTc nmc];
    spG.order=[kp kt km];
else
    spG.knots={Pknts,Tknts};
    spG.number=[nPc nTc];
    spG.order=[kp kt];
end
spG.dim=1;

if not(normflg)
   G_norm=1;
   Gp_norm=1;
   Gtt_norm=1;
end

% set up matrixes for solution of coef=A\B
% a number of different setups depending on all the choices for "data" and regularization
if reg_flg
    if m_flg
        if (rhoflg && Cpflg)
            if Gvflg          
                A=[G_norm*RdatG(id_rb,:);weight(1)*Gp_norm*RdatV(id_incl,:); weight(2)*Gtt_norm*RdatC(id_incl,:); norm_facP*lam(1)*RregP;norm_facT*lam(2)*RregT;norm_facm*lam(3)*Rregm];
                B=[G_norm*G(id_rb(:));weight(1)*Gp_norm*Gp(id_incl(:));weight(2)*Gtt_norm*Gtt(id_incl(:));zeros(3*nr,1)];
            else
                A=[G_norm*RdatG(id_incl,:);weight(1)*Gp_norm*RdatV(id_incl,:); weight(2)*Gtt_norm*RdatC(id_incl,:); norm_facP*lam(1)*RregP;norm_facT*lam(2)*RregT;norm_facm*lam(3)*Rregm];
                B=[G_norm*G(id_incl);weight(1)*Gp_norm*Gp(id_incl(:));weight(2)*Gtt_norm*Gtt(id_incl(:));zeros(3*nr,1)];
            end
        elseif (rhoflg && not(Cpflg))
            A=[G_norm*RdatG(id_incl,:);weight(1)*Gp_norm*RdatV(id_incl,:); norm_facP*lam(1)*RregP;norm_facT*lam(2)*RregT;norm_facm*lam(3)*Rregm];
            B=[G_norm*G(id_incl);weight(1)*Gp_norm*Gp(id_incl(:));zeros(2*nr,1)];      
        elseif (not(rhoflg) && Cpflg)
            A=[G_norm*RdatG(id_incl,:); weight(2)*Gtt_norm*RdatC(id_incl,:); norm_facP*lam(1)*RregP;norm_facT*lam(2)*RregT;norm_facm*lam(3)*Rregm];
            B=[G_norm*G(id_incl);weight(2)*Gtt_norm*Gtt(id_incl(:));zeros(3*nr,1)];        
        elseif (not(rhoflg) && not(Cpflg))
            A=[G_norm*RdatG(id_incl,:);norm_facP*lam(1)*RregP;norm_facT*lam(2)*RregT;norm_facm*lam(3)*Rregm];
            B=[G_norm*G(id_incl);zeros(3*nr,1)];      
        end        
    else
        if (rhoflg && Cpflg)
            if Gvflg          
                A=[G_norm*RdatG(id_rb,:);weight(1)*Gp_norm*RdatV(id_incl,:); weight(2)*Gtt_norm*RdatC(id_incl,:); norm_facP*lam(1)*RregP;norm_facT*lam(2)*RregT];
                B=[G_norm*G(id_rb(:))';weight(1)*Gp_norm*Gp(id_incl(:));weight(2)*Gtt_norm*Gtt(id_incl(:));zeros(2*nr,1)];
            else
                A=[G_norm*RdatG(id_incl,:);weight(1)*Gp_norm*RdatV(id_incl,:); weight(2)*Gtt_norm*RdatC(id_incl,:); norm_facP*lam(1)*RregP;norm_facT*lam(2)*RregT];
                B=[G_norm*G(id_incl);weight(1)*Gp_norm*Gp(id_incl(:));weight(2)*Gtt_norm*Gtt(id_incl(:));zeros(2*nr,1)];
            end
        elseif (rhoflg && not(Cpflg))
            A=[G_norm*RdatG(id_incl,:);weight(1)*Gp_norm*RdatV(id_incl,:); norm_facP*lam(1)*RregP;norm_facT*lam(2)*RregT];
            B=[G_norm*G(id_incl);weight(1)*Gp_norm*Gp(id_incl(:));zeros(2*nr,1)];      
        elseif (not(rhoflg) && Cpflg)
            A=[G_norm*RdatG(id_incl,:); weight(2)*Gtt_norm*RdatC(id_incl,:); norm_facP*lam(1)*RregP;norm_facT*lam(2)*RregT];
            B=[G_norm*G(id_incl);weight(2)*Gtt_norm*Gtt(id_incl(:));zeros(2*nr,1)];        
        elseif (not(rhoflg) && not(Cpflg))
            A=[G_norm*RdatG(id_incl,:);norm_facP*lam(1)*RregP;norm_facT*lam(2)*RregT];
            B=[G_norm*G(id_incl);zeros(2*nr,1)];      
        end
    end
else
    if (rhoflg && Cpflg)
        if Gvflg          
            A=[G_norm*RdatG(id_rb,:);weight(1)*Gp_norm*RdatV(id_incl,:); weight(2)*Gtt_norm*RdatC(id_incl,:)];
            B=[G_norm*G(id_rb(:));weight(1)*Gp_norm*Gp(id_incl(:));weight(2)*Gtt_norm*Gtt(id_incl(:))];
        else
            A=[G_norm*RdatG(id_incl,:);weight(1)*Gp_norm*RdatV(id_incl,:); weight(2)*Gtt_norm*RdatC(id_incl,:)];
            B=[G_norm*G(id_incl);weight(1)*Gp_norm*Gp(id_incl(:));weight(2)*Gtt_norm*Gtt(id_incl(:))];
        end        
    elseif (rhoflg && not(Cpflg))
        if Gvflg,error(' Cp is needed in this case'),end
        A=[G_norm*RdatG(id_incl,:);weight(1)*Gp_norm*RdatV(id_incl,:)];
        B=[G_norm*G(id_incl);weight(1)*Gp_norm*Gp(id_incl(:))];        
    elseif (not(rhoflg) && Cpflg)
        if Gvflg,error(' rho is needed in this case'),end
        A=[G_norm*RdatG(id_incl,:);weight(2)*Gtt_norm*RdatC(id_incl,:)];
        B=[G_norm*G(id_incl);weight(2)*Gtt_norm*Gtt(id_incl(:))];
    elseif (not(rhoflg) && not(Cpflg))
        if Gvflg,error(' rho and Cp are needed in this case'),end
        A=G_norm*RdatG(id_incl,:);
        B=G_norm*G(id_incl);
    end
end

switch algorthm
    case 'LeastSqr'
        if ndim==3
           spG.coefs=reshape(((A'*A)\(A'*B))',nPc,nTc,nmc);
        else
           spG.coefs=reshape(((A'*A)\(A'*B))',nPc,nTc);
        end
    case 'NormlEqn'
        if ndim==3
           spG.coefs=reshape((A\B),nPc,nTc,nmc); % solve normal equations
        else
           spG.coefs=reshape((A\B),nPc,nTc); % solve normal equations
        end
end

% alternative solutions: (try both compare - typically does not change the results)
%spG.coefs=reshape(((A'*A)\(A'*B))',nPc,nTc);  % least squares solution

spG.knots{1}=spG.knots{1}/1e6;  % convert back to MPa

% end of G spline creation

function D=makeCmn(A,B,C)
% This function combines collocation matrixes A(i,j) B(k,l) C(p,q) for 
% three dimensions (associated with a 4D hypersurface tensor spline)
% into a single matrix of size i*k*p by j*l*q. The idea is that the 
% tensor spline calculation is given as:
%    sum sum sum (ABC*a)  where A,B, and C are collocation matrix  and a is
%    an array of coefficients. 
% This can be remapped into a single matrix multiply D*a(:)
%    where D is the matrix returned here.
%  If working in 3D pass only matrixes A and B to this function
%  If A and B are row vectors, makeCmn returns a
%  single row corresponding to a single "data" point.
%
% This is mostly a problem of keeping track of how matrixes are indexed. 
% Below is the algorithm using full matrixes and a full set of "for" loops
% [ai,aj]=size(A);
% [bk,bl]=size(B);
% [cp,cq]=size(C);
% for p=1:cp
%   for k=1:bk      
%     for i=1:ai
%          mm= i + (k-1)*ai + bk*ai*(p-1);
%          a=A(i,:);  %one loop over an index can be vectorized
%          for q=1:cq  
%            for l=1:bl
%                 nnc=(l-1)*aj + bl*aj*(q-1);
%                 D(mm,nnc+(1:aj))=a*B(k,l)*C(p,q);
%            end
%          end
%      end
%   end
% end
% % % 
% actual calculation is shortened since the matrixes are very sparse
% JMB 2011

% Below I take advantage of sparseness - 
[ai,aj]=size(A);
[bk,bl]=size(B);
%Determine whether there is another dimension to the problem
switch nargin
    case 3
     [cp,cq]=size(C);
    case 2
     C=1;
     cp=1;cq=1;
end

% determine which elements are non-zero since A,B, and C are sparse
[ii,jj,sa]=find(A);
[kk,ll,sb]=find(B);
[pp,qq,sc]=find(C);

%Determine the sizes of the matrixes
nR= ai*bk*cp;  
nC= aj*bl*cq;
ni=length(ii);
nk=length(kk);
np=length(pp);
nmm=ni*nk*np;

% preallocate vectors
mm=zeros(nmm,1);
nn=zeros(nmm,1);
v=zeros(nmm,1);

%Do the work
for p=1:np
 for k=1:nk
     im= (kk(k)-1)*ai + bk*ai*(pp(p)-1);
     in= (ll(k)-1)*aj + bl*aj*(qq(p)-1);
     count=(k-1)*ni +ni*nk*(p-1);
     for i=1:ni 
       mm(count+i)= ii(i) + im;  
       nn(count+i)= jj(i) + in;
       v(count+i)=sa(i)*sb(k)*sc(p);
     end
 end
end

%Assemble the resultant sparse matrix
D=sparse(mm,nn,v,nR,nC);


function wt=area_wt(PT)
% Usage: wt=area_wt(PT) where PT is a cell containting the P and T grid points
% this subfunction assigns a weight to every regularization point based on
% the square root of the volume encompassed by that point
%  ie  sum of wt^2 is 1
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
 
function xg=mkgrid(x,nReg)
% this subfunction puts nreg regularization points between every control point
xg=x(:)';
if nReg>1
    dx=diff(xg);
    nx=length(xg);
    nt=(nx-1)*nReg+1;
    xtmp=zeros(1,nt);
    xtmp(1:nx)=xg;
    nstrt=nx+1;
    for i=2:nReg
        tmp=x(1:end-1)+(i-1)*dx/nReg;
        xtmp(nstrt:(nstrt+nx-2))=tmp;
        nstrt=nstrt+nx-1;
    end
    xg=sort(xtmp);
end

