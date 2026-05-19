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
%    kntflg default to 'c' for use of control points, 'k' to input the knots directly
%    lam is a vector of smoothing (damping) parameters (one for each independent variable) 
%           if lam is not provided, the values are set to 1
%    RegFac is a vector of integers for the regularization oversampling (number of regularization pts between each control point)
%           if RegFac is not included in options, the default  is 1 for all dimensions     
%           if RegFac is included but empty, do a least square fit without regularization
%                (likely to fail unless data coverage is adequate)      
%    ordr is a vector of integers defining the spline order for each independent variable
%          if ordr is not provided all are set to 4  (cubic splines)
%    mdrv is a vector of integers defining which derivative is to be  minimized
%          if mdrv is not provided set all to 2
%    mask is a list of ones and nan's to specify which data to exclude from the fit
%        if mask is not provided, all values of mask are set to 1
%    algorithm  if present, use Least Square fitting (default is Normal equations)
%
% This function requires the following external files:
%   collocate.m, makeCmn.m, sp_val.m, setknt.m
%   The subfunctions area_wt and mkgrid remain nested here.
%         JMB 2015-2026





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
    if (npts==1 || nXvar==1)
        X=X(:);
        nXvar=1;
        flg_grd=1;
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
         Xc{i}=X{i};   % X is a cell when flg_grd=1; was X(:,i) which is invalid
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
        flgReg=0;      % no regularization -- least-squares fit only
        RegFac=[];
    else
       RegFac=input_options.RegFac;
       flgReg=1;
    end
else
   flgReg=1;           % default: regularize with RegFac=1 in each dimension
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


tmp=std(Y(n_incld));
if(isnan(tmp)|| not(isfinite(tmp)) ),error('One or more input values to be fit are not finite'),end 
mdrv_indx=mdrv+1; %increment derivative by one to allign with indexing of basis functions from collocate

if (nC~=nXvar) 
    error('numer of Control point variables is not the same as the number of independent variables'),
end

if nC>3, error('numer of Control point variables is greater than 3 - not supported'),end

if nXvar==1
    xmin=min(X);
    xmax=max(X);
    if(Xc(1)>xmin || Xc(end)<xmax) 
        txt=sprintf('xmin=%4.1g xmax=%4.1g first knot=%4.1g last knot=%4.1g',xmin,xmax,Xc(1),Xc(end));
        error(['Data lie beyond the range of knots ' txt])
    end
else
    if flg_grd
        for i=1:nXvar
            xmin=min(X{i});
            xmax=max(X{i});
            if(Xc{i}(1)>xmin || Xc{i}(end)<xmax) 
                txt=sprintf('xmin=%4.1g xmax=%4.1g first knot=%4.1g last knot=%4.1g',xmin,xmax,Xc{i}(1),Xc{i}(end));
                error(['Data lie beyond the range of knots ' txt])
            end
        end
    else
        for i=1:nXvar 
              xmin=min(X(:,i));
              xmax=max(X(:,i));
              if(Xc{i}(1)>xmin || Xc{i}(end)<xmax) 
                txt=sprintf('xmin= %4.1g xmax= %4.1g first knot= %4.1g last knot= %4.1g',xmin,xmax,Xc{i}(1),Xc{i}(end));
                error(['Data lie beyond the range of knots ' txt])
              end         
        end
       
    end
end

kntflg='c';
if isfield(input_options,'kntflg')
    kntflg=input_options.kntflg;
end

nspl=1;
nnumber=zeros(1,nC);
switch kntflg
    case 'c'
        if(iscell(Xc))
            for i=1:nC
            nnumber(i)=length(Xc{i});
            nspl=nspl*nnumber(i);
            end
        else
           nnumber=length(Xc);
           nspl=nnumber;
        end
    case 'k'
        if(iscell(Xc))
            for i=1:nC           
            nnumber(i)=length(Xc{i})-2+ordr(i);
            nspl=nspl*nnumber(i);
            end
        else
           nnumber=length(Xc)-2+ordr(i);
           nspl=nnumber;
        end
end



if(isscalar(Uncert))
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
Wdm=norm(Wd,2);  % root mean square uncertainty
if(~isfinite(Wdm)),error('A data uncertainty is zero - needs to be fixed'), end

% if flg_grd
%    if (nY~=nXvar & nC~=1)
%        error('Y matix does not have the correct dimensions');
%    end
% end

switch kntflg
    case 'c'
       switch nC
         case 1
           Xknts = setknt(Xc(:)',ordr(1)) ;
         case 2
            Xknts = setknt(Xc{1}(:)',ordr(1)) ;
            Yknts = setknt(Xc{2}(:)',ordr(2)) ;
         case 3
            Xknts = setknt(Xc{1}(:)',ordr(1)) ;
            Yknts = setknt(Xc{2}(:)',ordr(2)) ;
            Zknts = setknt(Xc{3}(:)',ordr(3)) ;
       end
    case 'k'      
        switch nC
             case 1
               Xknts = [repmat(Xc{1}(1),1,ordr(1)) Xc{1}(2:end-1) repmat(Xc{1}(end),1,ordr(1))] ;
             case 2
                 Xknts = [repmat(Xc{1}(1),1,ordr(1)) Xc{1}(2:end-1) repmat(Xc{1}(end),1,ordr(1))] ;
                 Yknts = [repmat(Xc{2}(1),1,ordr(2)) Xc{2}(2:end-1) repmat(Xc{2}(end),1,ordr(2))] ;
             case 3
                 Xknts = [repmat(Xc{1}(1),1,ordr(1)) Xc{1}(2:end-1) repmat(Xc{1}(end),1,ordr(1))] ;
                 Yknts = [repmat(Xc{2}(1),1,ordr(2)) Xc{2}(2:end-1) repmat(Xc{2}(end),1,ordr(2))] ;
                 Zknts = [repmat(Xc{3}(1),1,ordr(3)) Xc{3}(2:end-1) repmat(Xc{3}(end),1,ordr(3))] ;
        end
end

% Creat the data co-location matrixes - different cases for different setups - 1
% 2 or 3 independent variables, gridded or scattered data, with or without
% parallel computing

if flg_grd
    switch nC
        case 1
            Rdat=collocate(Xknts,ordr(1),X,0);
        case 2
            Rdat=makeCmn(collocate(Xknts,ordr(1),X{1},0),collocate(Yknts,ordr(2),X{2},0));
        case 3
            Rdat=makeCmn(collocate(Xknts,ordr(1),X{1},0),collocate(Yknts,ordr(2),X{2},0),collocate(Zknts,ordr(3),X{3},0));
    end
else
    % Scattered data: collocate is vectorised over all points simultaneously,
    % so we evaluate all sites at once instead of looping point-by-point.
    % collocate requires sorted input, so we sort each dimension independently,
    % evaluate, then restore the original point ordering.
    switch nC
        case 1
            [x_srt, idx_x]   = sort(X(:,1));
            [~,     inv_x]   = sort(idx_x);
            col_srt          = collocate(Xknts, ordr(1), x_srt, 0);
            Rdat             = col_srt(inv_x, :, 1);   % (npts x nCoefs)

        case 2
            tmpx = X(:,1);  tmpy = X(:,2);
            [x_srt, idx_x] = sort(tmpx);  [~, inv_x] = sort(idx_x);
            [y_srt, idx_y] = sort(tmpy);  [~, inv_y] = sort(idx_y);
            colX_srt = collocate(Xknts, ordr(1), x_srt, 0);   % (npts x nX x 1)
            colY_srt = collocate(Yknts, ordr(2), y_srt, 0);   % (npts x nY x 1)
            % Restore original point ordering before assembling per-point rows
            colX_all = colX_srt(inv_x, :, 1);   % (npts x nX)
            colY_all = colY_srt(inv_y, :, 1);   % (npts x nY)
            Rdat = sparse(npts, nspl);
            for k = 1:npts
                Rdat(k,:) = makeCmn(colX_all(k,:), colY_all(k,:));
            end

        case 3
            tmpx = X(:,1);  tmpy = X(:,2);  tmpz = X(:,3);
            [x_srt, idx_x] = sort(tmpx);  [~, inv_x] = sort(idx_x);
            [y_srt, idx_y] = sort(tmpy);  [~, inv_y] = sort(idx_y);
            [z_srt, idx_z] = sort(tmpz);  [~, inv_z] = sort(idx_z);
            colX_srt = collocate(Xknts, ordr(1), x_srt, 0);
            colY_srt = collocate(Yknts, ordr(2), y_srt, 0);
            colZ_srt = collocate(Zknts, ordr(3), z_srt, 0);
            colX_all = colX_srt(inv_x, :, 1);
            colY_all = colY_srt(inv_y, :, 1);
            colZ_all = colZ_srt(inv_z, :, 1);
            Rdat = sparse(npts, nspl);
            for k = 1:npts
                Rdat(k,:) = makeCmn(colX_all(k,:), colY_all(k,:), colZ_all(k,:));
            end
    end
end

Nd=1/norm(Rdat,'Fro');

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
        colX=collocate(Xknts,ordr(1),Xg,mdrv_indx(1)-1);
        X2X=colX(:,:,mdrv_indx(1));
        [nr,~]=size(X2X);
        NX=1/norm(X2X,'Fro');
      case 2
          Xg=mkgrid(Xc{1},RegFac(1));
          Yg=mkgrid(Xc{2},RegFac(2));
          colX=collocate(Xknts,ordr(1),Xg,mdrv_indx(1)-1);
          colY=collocate(Yknts,ordr(2),Yg,mdrv_indx(2)-1);
          wtreg=area_wt({Xg,Yg});
          X2X=makeCmn(colX(:,:,mdrv_indx(1)),colY(:,:,1));
          Y2Y=makeCmn(colX(:,:,1),colY(:,:,mdrv_indx(2)));
         [nr,~]=size(X2X);
          nr=2*nr;
          NX=1/norm(X2X,'Fro');
          NY=1/norm(Y2Y,'Fro');
      case 3

        Xg=mkgrid(Xc{1},RegFac(1));
        Yg=mkgrid(Xc{2},RegFac(2));
        Zg=mkgrid(Xc{3},RegFac(3));
        wtreg=area_wt({Xg,Yg,Zg});
        colX=collocate(Xknts,ordr(1),Xg,mdrv_indx(1)-1);
        colY=collocate(Yknts,ordr(2),Yg,mdrv_indx(2)-1);
        colZ=collocate(Zknts,ordr(3),Zg,mdrv_indx(3)-1);
        X2X=makeCmn(colX(:,:,mdrv_indx(1)),colY(:,:,1),colZ(:,:,1));
        Y2Y=makeCmn(colX(:,:,1),colY(:,:,mdrv_indx(2)),colZ(:,:,1));
        Z2Z=makeCmn(colX(:,:,1),colY(:,:,1),colZ(:,:,mdrv_indx(3)));
        [nr,~]=size(X2X);
        nr=3*nr;
        NX=1/norm(X2X,'Fro');
        NY=1/norm(Y2Y,'Fro');
        NZ=1/norm(Z2Z,'Fro');
    end
end



% set up the problem to solve: coef=A\B  lots of different cases
if flgReg
    switch nC
        case 1
            if maskflg
                A=[Nd*Rdat(n_incld,:).*Wd(:);NX*lam(1)*(wtreg(:)*ones(1,nspl)).*X2X];
                Y=Y(:);
                B=[Nd*Y(n_incld(:)).*Wd(:);zeros(nr,1)];

            else
                A=[Nd*Rdat.*Wd(:);NX*lam(1)*wtreg(:).*X2X];
                B=[Nd*Y(:).*Wd(:);zeros(nr,1)];
            end
            switch algorthm
                case 'LeastSqr'
                    sp.coefs = reshape((A'*A)\(A'*B), 1, nnumber);
                case 'NormlEqn'
                    sp.coefs = reshape(A\B, 1, nnumber);   % (1 x n) row: sp convention for d=1
            end
        case 2
            if maskflg

                A=[Nd*Rdat(n_incld,:).*Wd(:);NX*lam(1)*wtreg(:).*X2X;NY*lam(2)*wtreg(:).*Y2Y];
                B=[Nd*Y(n_incld(:)).*Wd(:);zeros(nr,1)];
            else

                A=[Nd*Rdat.*Wd(:);NX*lam(1)*wtreg(:).*X2X;NY*lam(2)*wtreg(:).*Y2Y];
                B=[Nd*Y(:).*Wd(:);zeros(nr,1)];              
            end

            switch algorthm
                case 'LeastSqr'
                    sp.coefs=reshape((A'*A)\(A'*B),nnumber(1),nnumber(2)); 
                case 'NormlEqn'
                    sp.coefs=reshape(A\B,nnumber(1),nnumber(2)); 
            end
         case 3
             if maskflg               
                A=[Nd*Rdat(n_incld,:).*Wd(:);NX*lam(1)*wtreg(:).*X2X;NY*lam(2)*wtreg(:).*Y2Y;NZ*lam(3)*wtreg(:).*Z2Z];
                B=[Nd*Y(n_incld).*Wd(:);zeros(nr,1)];
             else
                A=[Nd*Rdat.*Wd(:);NX*lam(1)*wtreg(:).*X2X;NY*lam(2)*wtreg(:).*Y2Y;NZ*lam(3)*wtreg(:).*Z2Z];
                B=[Nd*Y(:).*Wd(:);zeros(nr,1)];
             end
             switch algorthm
                case 'LeastSqr'        
                     sp.coefs=reshape((A'*A)\(A'*B),nnumber(1),nnumber(2),nnumber(3));
                case 'NormlEqn'
                     sp.coefs=reshape(A\B,nnumber(1),nnumber(2),nnumber(3));
             end
    end
else
    switch nC
        case 1
            A=Nd*Rdat.*Wd(:);
            B=Nd*Y(:).*Wd(:);
            switch algorthm
                case 'LeastSqr'
                    sp.coefs = reshape((A'*A)\(A'*B), 1, nnumber);
                case 'NormlEqn'
                    sp.coefs = reshape(A\B, 1, nnumber);   % (1 x n) row: sp convention for d=1
            end
        case 2
            A=Nd*Rdat.*Wd(:);
            B=Nd*Y(:).*Wd(:);
            switch algorthm
                case 'LeastSqr'
                    sp.coefs=reshape(((A'*A)\(A'*B)),nnumber(1),nnumber(2));
                case 'NormlEqn'
                    sp.coefs=reshape(A\B,nnumber(1),nnumber(2));
            end
            
         case 3
            A=Nd*Rdat.*Wd(:);
            B=Nd*Y(:).*Wd(:);
            switch algorthm
                case 'LeastSqr'
                    sp.coefs=reshape(((A'*A)\(A'*B)),nnumber(1),nnumber(2),nnumber(3));
                case 'NormlEqn'
                    sp.coefs=reshape(A\B,nnumber(1),nnumber(2),nnumber(3));
            end
         
    end
end

% determine statististcs only for the range defined by the control points

%if flg_grd
 Yc=sp_val(sp,X);
% else
%  Yc=sp_val(sp,X);  % needs transpose because the convention here is columns of data and fnval requires rows of data
% end
% Ycc=sp_val(sp,Xc);
 
devs=Y(n_incld(:))-Yc(n_incld(:));
ppm=1e6*sqrt(sum(devs.^2./Y(n_incld(:)).^2)/ndat);

rms=sqrt(sum(devs.^2)/ndat);

chisqr=sum((devs(:)./Uncert(:)).^2)/ndat;

sp.Data.X=X;
sp.Data.Y=Y;
sp.Data.Yc=Yc;
%sp.Data.Ycontrl=Ycc;
sp.Data.Uncert=Uncert;
sp.Data.Xc=Xc;
sp.Data.lam=lam;
sp.Data.mdrv=mdrv;
sp.Data.RegFac=RegFac;
sp.Data.rms=rms;
sp.Data.devs=devs;
sp.Data.chisqr=chisqr;
sp.Data.ppm=ppm;
sp.revision=datetime;

% -----------------------------------------------------------------------
% makeCmn, collocate, and setknt are external standalone files.
% The subfunctions below (area_wt, mkgrid) are private to spdft.
% -----------------------------------------------------------------------

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
if nreg>1
for i=2:nreg
    xg=[xg xtmp(1:end-1)+(i-1)*dx/nreg];
end
end
xg=sort(xg);

