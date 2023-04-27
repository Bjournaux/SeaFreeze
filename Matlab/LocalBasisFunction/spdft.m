function sp=spdft(X,Y,dY,input_options)
% Fit b splines to z=f(x,y,...) data -ie a scalar value z as a function
%      of one, two, or three independent variables
%      with regularization on  each independent variable.  
%
%Usage:   sp=spdft(X,Y,dY,options)
% where: The input (X and Y) can be gridded or scattered and can be of variable
% (2, 3, or 4) dimensions.  Gridded data are entered as cells and scattered data as
% matrixes with columns for each independent variable.  
%        
%     X contains independent variables (matrix or cell), Y is the dependent variable
%      dY is the data uncertainty.  (1) If empty, use stdev(Y) as weight (2)
%      If scalar, assign same uncertainty to all data, (3) If same size as Y
%      individual uncertainties are assigned
%
% The following input options are entered in the structure "options"
%    Xc (vector or cell-if more than one independent variable)are the "control points" used to determine knot locations.
%             If Xc is not provided and X is a cell, Xc is set to X
%    lam is a vector of smoothing parameters (one for each independent variable) 
%           if lam is not provided, the values are set to 1
%    RegFac is a vector of integers for the regularization oversampling (number of regularization pts between each control point)
%           if RegFac is not included in options, the default  is 1 for all dimensions     
%           if RegFac is included but empty, do a least square fit without regularization
%                (likely to fail unless data coverage is adequate)      
%    ordr is a vector of integers defining the spline order for each independent variable
%          if ordr is not provided set all to 4
%    mdrv is a vector of integers defining which derivative is to be  minimized
%          if mdrv is not provided set all to 2
%    mask is a list of ones or nan to specify which data to exclude from the fit
%        if mask is not provided, set all values of mask to 1
%    algorithm  to use Least Square fitting (default is Normal equations)
%
% This function requires routines from the matlab spline (curve fitting) toolbox.  
%   Custom function mkCmn mkgrid area_wt are nested in this file.
%         JMB 2015

%test whether parallel computing is possible on this computer:
has_parfor = ~isempty(which('parfor'));


if nargin==3
    input_options=[];
end

Y=Y(:);  % force dependent variable to always be in a column

flg_grd=0;
if(iscell(X))
    flg_grd=1;
    [~,nXvar]=size(X);
    nX=zeros(1,nXvar);   
    for i=1:nXvar
        nX(i)=length(X{i});      
    end
else
    [npts,nXvar]=size(X);
    if npts==1
        X=X';
        nXvar=1;
    end
end

if isfield(input_options,'Xc')
    Xc=input_options.Xc;
    if(iscell(Xc))
      nC=length(Xc);
    else
      nC=1;
    end
else
    if ((flg_grd) && nXvar>1)
      for i=1:nXvar
         Xc{i}=X(:,i);
      end
    end
    nC=nXvar;
end

if isfield(input_options,'lam')
    lam=input_options.lam;
else
    lam=ones(1,nXvar);
end

if isfield(input_options,'RegFac')
    if isempty(input_options.RegFac)
        flgReg=0;
    else
       RegFac=input_options.RegFac;
       flgReg=1;
    end
else
   flgReg=1;  %just do least square fit data without regularization   
   RegFac=ones(1,nXvar);
end

if isfield(input_options,'ordr')
    ordr=input_options.ordr;
else
    ordr=4*ones(1,nXvar);
end

if isfield(input_options,'mdrv')
    mdrv=input_options.mdrv;
else
    mdrv=2*ones(1,nXvar);
end

if isfield(input_options,'mask')
    maskflg=1;
    mask=input_options.mask;
    n_incld=find(isfinite(mask(:)));
    ndat=length(n_incld);
    if(max(size(Y(:))~=size(mask(:)))~=0),error('mask size and size of input data are not compatible'),end
else
    maskflg=0;
    mask=ones(size(Y));
    ndat=length(Y(:));
    n_incld=1:ndat;
end

if isempty(dY)
    Uncert=std(Y(n_incld));   
else
    Uncert=dY(:);
end

if(isfield(input_options,'algorithm'))
    algorthm='LeastSqr';
else
    algorthm='NormlEqn';
end

NA=sqrt(ndat);   % weight for the total number of data being fit
tmp=std(Y(n_incld));
if(isnan(tmp)|| not(isfinite(tmp)) ),error('One or more input values to be fit are not finite'),end 
mdrv_ordr=mdrv+1; %increment derivative by one to allign with indexing of basis functions from spcol

if (nC~=nXvar) 
    error('numer of Control point variables is not the same as the number of independent variables'),
end

if nC>3, error('numer of Control point variables is greater than 3 - not supported'),end

if nXvar==1
    xmin=min(X);
    xmax=max(X);
    if(Xc(1)>xmin || Xc(end)<xmax), error('Data lie beyond the range of knots'),end
else
    if flg_grd
        for i=1:nXvar
            xmin=min(X{i});
            xmax=max(X{i});
            Xc{i}(1);
            Xc{i}(end) ;         
            if(Xc{i}(1)>xmin || Xc{i}(end)<xmax), error('Data lie beyond the range of knots'),end
        end
    else
        for i=1:nXvar
           xmin=min(X(:,i));
           xmax=max(X(:,i));

           if(Xc{i}(1)>xmin || Xc{i}(end)<xmax), error('Data lie beyond the range of knots'),end           
        end
    end
end

nspl=1;
nnumber=zeros(1,nC);
if(iscell(Xc))
    for i=1:nC
    nnumber(i)=length(Xc{i});
    nspl=nspl*nnumber(i);
    end
else
   nnumber=length(Xc);
   nspl=nnumber;
end

if(length(Uncert)==1)
    if maskflg
        Uncert=Uncert*ones(size(Y(n_incld)));
    else
        Uncert=Uncert*ones(ndat,1);
    end
else
    if maskflg
        Uncert=Uncert(n_incld);
    end
end

Wd=Uncert(:).^(-1);               % uncertainty weights
Wdm=norm(Wd,2)/sqrt(length(Wd));  % root mean square uncertainty

if flg_grd
   nY=length(size(Y));
   if nY~=nXvar
       error('Y matix does not have the correct dimensions');
   end
end

switch nC
     case 1
       Xknts = optknt(Xc(:)',ordr(1)) ;
     case 2
        Xknts = optknt(Xc{1}(:)',ordr(1)) ;
        Yknts = optknt(Xc{2}(:)',ordr(2)) ;
     case 3
        Xknts = optknt(Xc{1}(:)',ordr(1)) ;
        Yknts = optknt(Xc{2}(:)',ordr(2)) ;
        Zknts= optknt(Xc{3}(:)',ordr(3)) ;
end

% Creat the data co-location matrixes - different cases for different setups - 1
% 2 or 3 independent variables, gridded or scattered data, with or without
% parallel computing
if flg_grd
    switch nC
        case 1
            Rdat=spcol(Xknts,ordr(1),brk2knt(X,1));
        case 2
            Rdat=makeCmn(spcol(Xknts,ordr(1),brk2knt(X{1},1)),spcol(Yknts,ordr(2),brk2knt(X{2},1)));
        case 3
            Rdat=makeCmn(spcol(Xknts,ordr(1),brk2knt(X{1},1)),spcol(Yknts,ordr(2),brk2knt(X{2},1)),spcol(Zknts,ordr(3),brk2knt(X{3},1)));
    end
else    
    Rdat=sparse(npts,nspl);
    switch nC
       case 1
        if has_parfor 
         parfor k=1:npts
           Rdat(k,:)=spcol(Xknts,ordr(1),brk2knt(X(k,1),1));
         end
        else
         for k=1:npts
           Rdat(k,:)=spcol(Xknts,ordr(1),brk2knt(X(k,1),1));
         end
        end
       case 2
           tmpx=X(:,1);
           tmpy=X(:,2);
          if has_parfor 
           parfor k=1:npts
             Rdat(k,:)=makeCmn(spcol(Xknts,ordr(1),brk2knt(tmpx(k),1)),spcol(Yknts,ordr(2),brk2knt(tmpy(k),1)));
           end  
          else
           for k=1:npts
             Rdat(k,:)=makeCmn(spcol(Xknts,ordr(1),brk2knt(tmpx(k),1)),spcol(Yknts,ordr(2),brk2knt(tmpy(k),1)));
           end  
          end
          
       case 3
           tmpx=X(:,1);
           tmpy=X(:,2);
           tmpz=X(:,3);
           if has_parfor 
           parfor k=1:npts
             Rdat(k,:)=makeCmn(spcol(Xknts,ordr(1),brk2knt(tmpx(k),1)),spcol(Yknts,ordr(2),brk2knt(tmpy(k),1)),spcol(Zknts,ordr(3),brk2knt(tmpz(k),1)));
           end
           else
           for k=1:npts
             Rdat(k,:)=makeCmn(spcol(Xknts,ordr(1),brk2knt(tmpx(k),1)),spcol(Yknts,ordr(2),brk2knt(tmpy(k),1)),spcol(Zknts,ordr(3),brk2knt(tmpz(k),1)));
           end
           end
     end
end

   
sp.form='B-';
switch nC
    case 1
       sp.knots=Xknts;
       sp.number=nnumber;
       sp.order=ordr;
    case 2
       sp.knots={Xknts,Yknts};
       sp.number=nnumber;
       sp.order=ordr; 
    case 3
       sp.knots={Xknts,Yknts,Zknts};
       sp.number=nnumber;
       sp.order=ordr;
end
sp.dim=1;

% Generate collocation matrixes for regularization: different number of
% independent variables
if flgReg
    switch nC
      case 1
        Xg=mkgrid(Xc,RegFac(1));
        wtreg=area_wt(Xg);
        colX=spcol(Xknts,ordr(1),brk2knt(Xg,mdrv_ordr(1)));
        X2X=colX(mdrv_ordr(1):mdrv_ordr(1):end,:);
        [nr,~]=size(X2X);
      case 2
          Xg=mkgrid(Xc{1},RegFac(1));
          Yg=mkgrid(Xc{2},RegFac(2));
          colX=spcol(Xknts,ordr(1),brk2knt(Xg,mdrv_ordr(1)));
          colY=spcol(Yknts,ordr(2),brk2knt(Yg,mdrv_ordr(2)));
          wtreg=area_wt({Xg,Yg});
          X2X=makeCmn(colX(mdrv_ordr(1):mdrv_ordr(1):end,:),colY(1:mdrv_ordr(2):end,:));
          Y2Y=makeCmn(colX(1:mdrv_ordr(1):end,:),colY(mdrv_ordr(2):mdrv_ordr(2):end,:));
         [nr,~]=size(X2X);
          nr=2*nr;
      case 3
        Xg=mkgrid(Xc{1},RegFac(1));
        Yg=mkgrid(Xc{2},RegFac(2));
        Zg=mkgrid(Xc{3},RegFac(3));
        wtreg=area_wt({Xg,Yg,Zg});
        colX=spcol(Xknts,ordr(1),brk2knt(Xg,mdrv_ordr(1)));
        colY=spcol(Yknts,ordr(2),brk2knt(Yg,mdrv_ordr(2)));
        colZ=spcol(Zknts,ordr(3),brk2knt(Zg,mdrv_ordr(3)));
        X2X=makeCmn(colX(mdrv_ordr(1):mdrv_ordr(1):end,:),colY(1:mdrv_ordr(2):end,:),colZ(1:mdrv_ordr(3):end,:));
        Y2Y=makeCmn(colX(1:mdrv_ordr(1):end,:),colY(mdrv_ordr(2):mdrv_ordr(2):end,:),colZ(1:mdrv_ordr(3):end,:));
        Z2Z=makeCmn(colX(1:mdrv_ordr(1):end,:),colY(1:mdrv_ordr(2):end,:),colZ(mdrv_ordr(3):mdrv_ordr(3):end,:));
        [nr,~]=size(X2X);
        nr=3*nr;
    end
end

% set up the problem to solve: coef=A\B  lots of different cases
if flgReg
    switch nC
        case 1
            if maskflg
                A=[Rdat(n_incld,:).*(Wd(:)*ones(1,nspl));Wdm*NA*lam(1)*(wtreg(:)*ones(1,nspl)).*X2X]/Wdm;
                Y=Y(:);
                B=[Y(n_incld(:)).*Wd(:);zeros(nr,1)]/Wdm;
                sp.coefs=(A\B)';
            else
                A=[Rdat.*(Wd(:)*ones(1,nspl));Wdm*NA*lam(1)*(wtreg(:)*ones(1,nspl)).*X2X]/Wdm;
                B=[Y(:).*Wd(:);zeros(nr,1)]/Wdm;
                sp.coefs=(A\B); 
            end
        case 2
            if maskflg
                A=[Rdat(n_incld,:).*(Wd(:)*ones(1,nspl));Wdm*NA*lam(1)*(wtreg(:)*ones(1,nspl)).*X2X;Wdm*NA*lam(2)*(wtreg(:)*ones(1,nspl)).*Y2Y];
                B=[Y(n_incld(:)).*Wd(:);zeros(nr,1)];
                sp.coefs=reshape(A\B,nnumber(1),nnumber(2));
            else
                A=[Rdat.*(Wd(:)*ones(1,nspl));Wdm*NA*lam(1)*(wtreg(:)*ones(1,nspl)).*X2X;Wdm*NA*lam(2)*(wtreg(:)*ones(1,nspl)).*Y2Y];
                B=[Y(:).*Wd(:);zeros(nr,1)];
                sp.coefs=reshape(A\B,nnumber(1),nnumber(2)); 
            end
         case 3
             if maskflg               
                A=[Rdat(n_incld,:).*(Wd(:)*ones(1,nspl));Wdm*NA*lam(1)*(wtreg(:)*ones(1,nspl)).*X2X;Wdm*NA*lam(2)*(wtreg(:)*ones(1,nspl)).*Y2Y;Wdm*NA*lam(3)*(wtreg(:)*ones(1,nspl)).*Z2Z]/Wdm;
                B=[Y(n_incld).*Wd(:);zeros(nr,1)]/Wdm;
                sp.coefs=reshape(A\B,nnumber(1),nnumber(2),nnumber(3));
             else
                A=[Rdat.*(Wd(:)*ones(1,nspl));Wdm*NA*lam(1)*(wtreg(:)*ones(1,nspl)).*X2X;Wdm*NA*lam(2)*(wtreg(:)*ones(1,nspl)).*Y2Y;Wdm*NA*lam(3)*(wtreg(:)*ones(1,nspl)).*Z2Z]/Wdm;
                B=[Y(:).*Wd(:);zeros(nr,1)]/Wdm;
                sp.coefs=reshape(A\B,nnumber(1),nnumber(2),nnumber(3));
             end
    end
else
    switch nC
        case 1
            A=Rdat.*(Wd(:)*ones(1,nspl));
            B=Y(:).*Wd(:);
            switch algorthm
                case 'LeastSqr'
                    sp.coefs=((A'*A)\(A'*B));
                case 'NormlEqn'
                    sp.coefs=(A\B);
            end
        case 2
            A=Rdat.*(Wd(:)*ones(1,nspl));
            B=Y(:).*Wd(:);
            switch algorthm
                case 'LeastSqr'
                    sp.coefs=reshape(((A'*A)\(A'*B)),nnumber(1),nnumber(2));
                case 'NormlEqn'
                    sp.coefs=reshape(A\B,nnumber(1),nnumber(2));
            end
            
         case 3
            A=Rdat.*(Wd(:)*ones(1,nspl));
            B=Y(:).*Wd(:);
            switch algorthm
                case 'LeastSqr'
                    sp.coefs=reshape(((A'*A)\(A'*B)),nnumber(1),nnumber(2),nnumber(3));
                case 'NormlEqn'
                    sp.coefs=reshape(A\B,nnumber(1),nnumber(2),nnumber(3));
            end
         
    end
end

% determine statististcs only for the range defined by the control points

if flg_grd
 Yc=fnval(sp,X);
else
 Yc=fnval(sp,X')';  % needs transpose because the convention here is columns of data and fnval requires rows of data
end
 Ycc=fnval(sp,Xc);
 
devs=Yc(n_incld(:))-Y(n_incld(:));

rms=sqrt(sum(devs.^2)/ndat);

chisqr=sum((devs(:)./Uncert(:)).^2)/ndat;

sp.Data.X=X;
sp.Data.Y=Y;
sp.Data.Yc=Yc;
sp.Data.Ycontrl=Ycc;
sp.Data.Uncert=Uncert;
sp.Data.Xc=Xc;
sp.Data.lam=lam;
sp.Data.mdrv=mdrv;
sp.Data.RegFac=RegFac;
sp.Data.rms=rms;
sp.Data.devs=devs;
sp.Data.chisqr=chisqr;

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
% Below is the algorithm using a full set of for loops
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
%                 b=a*B(k,l)*C(p,q);
%                 D(mm,nnc+(1:aj))=b;
%            end
%          end
%      end
%   end
% end
% % % 
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

%Determine (1) which elements are non-zero (A,B,C are sparse) and (2) the sizes of the matrixes
nR= ai*bk*cp;  
nC= aj*bl*cq;
[ii,jj,sa]=find(A);
[kk,ll,sb]=find(B);
[pp,qq,sc]=find(C);
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
xtmp=x(:)';
xg=xtmp;
dx=diff(xg);
if nreg>1;
for i=2:nreg
    xg=[xg xtmp(1:end-1)+(i-1)*dx/nreg];
end
end
xg=sort(xg);
