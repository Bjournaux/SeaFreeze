function spG=spgft_P(Data,Options)
% create a local basis function representation (tensor b splines) for Gibbs energy - use G (optionally rho and Cp) as constraint and
% regularize in the non-physical region by requiring smooth mdrv'th derivatives.
% This solves the linear problem where G can be function of P and T or function of P T and composition (m)
%
%    Units for all quantities are MKS!! (J/kg, ,MPa, K kg/m^3, kg/mol)  
%      (internal factor of 10^6  converts input MPa to Pa for internal calculations and then converts back for output)
% 
%    G rho and Cp are grids of values with PTm containing cells with P, T, and composition values.  
%        Alternative is to enter G only along a base (reference pressure isobar)-
%         In that case G is a vector of values and Cp and rho are used to reconstruct G by colocation at higher pressures
%
%    The spline is on a rectangular (can be non-uniform) grid of points. In the real world, 
%    phase boundaries interfere and sections of the grid represent metastable regions where 
%    no reasonable thermodynamic values are available.  Thus, to solve the problem on the grid (insure that the
%    matrixes are adequately conditions and not too rank-deficit) regularization is applied 
%        - assert that derivatives of the G surface are adequately "smooth" where data are not available. 
%
% function call:
%      sp_G=spgft(Data,Options)
%   where Data is stucture with
%    input "MW" the molecular weight for a compositional fit
%    input "PTm" a cell containing P and T or P T and m that defines the grid
%    input "G" is either gridded data at PTm values or Tm values for P=Pr
%    input "P_ref" is the pressure applying to vector G as a function of T (default P_ref=0.1 MPa)
%    inputs "rho" and "Cp" are either gridded data at PTm values or Tm values for P=Pr
%    input "mask" has ones where input thermodynamic properties are defined and NaN where they are not (Default is all ones)
%
%    all remaining parameters are in a structure called "Options" and have predefined default values 
%       "PTmc" is a cell containing the control points for the spline (Default is to use the PTm grid)              
%       "lam" is vector of smoothing parameters (default is all ones )
%       "weight" is vector of weights for the relative importance of densities and Cp relative to G (default is [1 1]);
%       "mdrv" is vector of the derivatives to use in smoothing (default is
%       [4 4 3] for smooth derivatives associated with compressibilites, specific heats, and chemicals potentials)
%       "ordr" is vector of integer order's for the spline  (default is [6 6 4])
%       "nReg" sets the number of regularization points between each control point (default is  1)
%       "normflg" =1 to rescale data based on their standard deviations (default is no scaling)
%       "algorithm" =1  to use Least Square fitting (default is Normal equations)
%
%  Spline toolbox functions called:
%           aptknt
%           spcol
%           brk2knt
%  custom (embeded) functions:
%              mkCmn:    maps multidimensional sparse array into matrix
%              wt_area:  determines relative area of points on non-uniform grid
%           and mkgrid:  determines locations of regularization points
%  Work remains to allow normalization schemes to better compensate for a
%  wide range of problems.  The current version still requires large
%  changes in damping and weights for different "size" problems - the goal
%  is to find a system that allows weights and damping near one for any
%  problem.
%  JMB 2015-2018 - 

if(isfield(Data,'PTm'))
    PTm=Data.PTm;
else
    error('Input "Data" needs a "PTm" field to define pressure-temperature-(optional)concentration grid')
end
ndim=length(PTm);
m_flg=0;
if ndim==3
    m_flg=1;
    if(isfield(Data,'nu'))
       SpG.MW=Data.nu;
    else
        error('Input "Data" needs a "nu" field for mixtures (count of species in solution) ')
    end
end

if(isfield(Data,'MW'))
   SpG.MW=Data.MW;
else
        error('Input "Data" needs a "MW" field (the molecular weights-kg/mol - scalar for single component and vector for mixtures) ')
end



if(isfield(Data,'rho'))
    rho=Data.rho;
else
    error('Input "Data" needs a "rho" field to define density grid')
end

if(isfield(Data,'Cp'))
    Cp=Data.Cp;
else
    error('Input "Data" needs a "Cp" field to define specific heat grid')
end

if(isfield(Data,'G'))
    G=Data.G;
else
    error('Input "Data" needs a "G" field to define Gibbs energy grid')
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


if Gvflg
  if(isfield(Data,'P_ref'))
    P_ref=Data.P_ref;
  else
    P_ref=0.1;
  end
else
 if(isfield(Data,'P_ref'))
     error('P_ref incorrectly set for G as a surface in P and T')
 end
end

if(isfield(Data,'mask'))
    mask=Data.mask;
else
    mask=ones(size(rho));
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
P_ref=P_ref*1e6;
T=PTm{2};
T=T(:)';
if m_flg
    m=PTm{3};
    m=m(:)';
end

if(isfield(Options,'PTmc'))
    PTmc=Options.PTmc;
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

if(isfield(Options,'normflg'))
    normflg=1;
else
    normflg=0;
end



reg_flg=1;
if(isfield(Options,'nReg'))
    nReg=Options.nReg;
    if(isempty(nReg))
      reg_flg=0;
    end
else
    nReg=ones(1,ndim);
end

if(isfield(Options,'ordr'))
    ordr=Options.ordr;
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

if(isfield(Options,'mdrv'))
    mdrv=Options.mdrv;
else
    if ndim==2
      mdrv=4*ones(1,2);
    elseif ndim==3
        mdrv=[4 4 3];
    end
end

if(isfield(Options,'lam'))
    lam=Options.lam;
else
    if ndim==2
      lam=[1 1];
    elseif ndim==3
        lam=[1 1 1];
    end
end


if(isfield(Options,'weight'))
    weight=Options.weight;
else
    weight=[1 1];
end

if(isfield(Options,'algorithm'))
    algorthm='LeastSqr';
else
    algorthm='NormlEqn';
end

% determine number of data and control points
nP=length(P);

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


m_rb=mask(1,:,:);
id_rb=find(not(isnan(m_rb(:))));

n_dat=length(id_incl);
norm_fac=sqrt(n_dat);  % will weight data vs regularization by square root of number of data points


if Gvflg % handle G only at reference pressure or G for all P and T and m
    tmp=std(G(id_rb));
    if tmp>0
     G_norm=tmp^-1; % use standard deviation of G at 1 bar for normalization
    else
        G_norm=1;
    end
else
    G_norm=std(G(id_incl(:)))^-1; % use standard deviation of G(P,T) for normalization
end

% rho and Cp expressed as linear derivatives of G and determine standard
% deviations of the "data" derivatives
if rhoflg 
    Gp=rho.^-1;
    Gp_norm=std(Gp(id_incl(:)))^-1; % use standard deviation of Gp for normalization
end    
if Cpflg
    if m_flg
        Gtt=-Cp./repmat(T(:)',nP,1,nm);
        Gtt_norm=std(Gtt(id_incl(:)))^-1;
    else
        id=find(T==0); 
          Gtt=-Cp./repmat(T(:)',nP,1);
          if(not(isempty(id)))
              Gtt(id)=0;
          end
          Gtt_norm=std(Gtt(id_incl(:)))^-1 ;
            
    end   % use standard deviation of Gtt for normalization
end

% check that all input data are finite
%   check  G_norm, Gp_norm,and Gtt_norm which should be finite if the input matrices are
if (not(isfinite(G_norm))), error('Some input G values are not finite'), end
if rhoflg
 if (not(isfinite(Gp_norm))), error('Some input rho values are not finite'), end
end
if Cpflg
  if (not(isfinite(Gtt_norm))), error('Some input Cp values are not finite'), end
end

% check the input normflg - set normalizations to 1 if not normalizing
if not(normflg)
   G_norm=1;
   Gp_norm=1;
   Gtt_norm=1;
end

% determine knot locations
Pknts=aptknt(Pc(:)',kp);
Tknts=aptknt(Tc(:)',kt);

Tcol = spcol(Tknts,kt,brk2knt(T,kt));
Pcol = spcol(Pknts,kp,brk2knt(P,kp));
Pcol_ref=spcol(Pknts,kp,brk2knt(P_ref,1));


if m_flg
    mknts=aptknt(mc(:)',km);
    mcol = spcol(mknts,km,brk2knt(m,km));
end

% create the matrixes of basis functions
% depends on whether G is defined everywhere or just at the reference pressure

if Gvflg
    if m_flg
        RdatG=makeCmn(Pcol_ref,Tcol(1:kt:end,:),mcol(1:km:end,:));
    else
        RdatG=makeCmn(Pcol_ref,Tcol(1:kt:end,:));
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
% then weight by the volume contained around each regularization point
% (accounting for possibly non-uniform grids)
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
% %NV=norm(RdatV,1);
% NC=norm(RdatC,1);
% NRP=norm(RregP,1);
% NRT=norm(RregT,1);
% % NRm=norm(Rregm,1);
% Gp_norm=Gp_norm*NG/NV;
% Gtt_norm=Gtt_norm*NG/NC;
% norm_facP=NG/NRP;
% norm_facT=NG/NRT;
% %norm_facm=NG/NRm;

% don't scale the regularization matrixes
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

% set up matrixes for the solution that gives the coefficients as A\B (normal equations) or (A'*A)\(A'*B) (least square solution)
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
            B=[G_norm*G(id_rb(:))';weight(1)*Gp_norm*Gp(id_incl(:));weight(2)*Gtt_norm*Gtt(id_incl(:))];
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

spG.knots{1}=spG.knots{1}/1e6;  % convert pressure back to MPa

% for use elsewhere, keep a copy of the control points used to make this
% spline.
if m_flg
 spG.PTmc={Pc,Tc,mc};
else
  spG.PTmc={Pc,Tc};
end  
spG.version=datetime;
% end of G spline creation - start embeded functions

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
% Usage: wt=area_wt(PT) where PT is a cell containting the P and T (+/- m) grid points
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

