function out=sp_val(sp,derv,x)
% function sp_val is offered to allow use of LBF representations without
% installation of the optional Curevefitting Toolbox (that includes spline functions).  
% sp_val substitutes for the combinations of fnval and fnder plus various utility 
% functions in the spline toolbox.  sp_val is "self contained" and calls no other
% specialized functions. The algorithms are those given by de Boor.  
%
%  Usage:   out = sp_val(LBF,derv,x) or = sp_val(LBF,x) if derivatives are not needed
%
%  where:    "out" is a vector or matrix of values (depending on the form of"x"
%            "LBF" is a b spline (tensor or univarient)
%            "derv" is vector of requested derivatives (one element per degree of freedom)
%            "x" is either a cell containing the independent variables (gridded output)
%                 or an array of scattered data.
%             sp_val(LBF,[0 0 ..],x) is equivalent to sp_val(LBF,x)
%  example: for an LBF in pressure (0-100 MPa), temperature (250-350 K) and concentration (0-1 M),
%           PTm={0:10:100,250:10:350,0:.1:1} and for first derivatives of
%           all independent variables:
%                  out=sp_val(LBF,[1 1 1],PTm);
%
% no error checking is included - user must understand how to use this function. 
%
% a major limitation of this function is that scattered data are calculated
% in a loop using the b-form splines (parfor speeds it up).  For scattered data, 
% the Mathworks provided fnval convert the b-splines to a pp-form and then does the 
% calculation using the pp-form. To do this requires more implementation of various
% spline manipulations. The slower method used here is about 3-5 times slower.
% 7/2019



if nargin==2,x=derv;derv=zeros(size(x(1,:)));end
   sp=fndr(sp,derv);  
   out = spvl(sp,x);
end
       
function out = spvl(sp,x)
% based on SPVAL that evaluates a function in B-form.

if iscell(sp.knots)  % we are dealing with a tensor product spline
   [t,a,n,~,d] = spbk(sp); m = length(t);
   nd=length(x);
   if iscell(x)  % evaluation on a mesh
     v = a; sizev = [d,n]; nsizev = zeros(1,m);
      for i=m:-1:1
         nsizev(i) = length(x{i}(:));
         v = reshape(...
         spv(spmk(t{i},reshape(v,prod(sizev(1:m)),sizev(m+1))), ...
                 x{i}),   [sizev(1:m),nsizev(i)]);
         sizev(m+1) = nsizev(i);
         if m>1
            v = permute(v,[1,m+1,2:m]); sizev(2:m+1) = sizev([m+1,2:m]);
         end
      end
      if d>1
         out = reshape(v,[d,nsizev]);
      else
         out = reshape(v,nsizev);
      end
 
   else          % evaluation at scattered points;
                 % this will eventually be done directly here.
      out=zeros(nd,1);
      parfor ii=1:nd
         sizev = [d,n]; nsizev = ones(1,m);
         v = a;
      for i=m:-1:1
         v = reshape(...
         spv(spmk(t{i},reshape(v,prod(sizev(1:m)),sizev(m+1))), ...
                 x(ii,i)),   [sizev(1:m),1]);
         sizev(m+1) = 1;
         if m>1
            v = permute(v,[1,m+1,2:m]); sizev(2:m+1) = sizev([m+1,2:m]);
         end
      end
      if d>1
         v = reshape(v,[d,nsizev]);
      else
         v = reshape(v,nsizev);
      end
      out(ii)=v;
      end
                                   
   end
else                 % we are dealing with a univariate spline 
   v = spv(sp,x);
end
end

function v = spv(sp,x)
% based on SPVAL1 that evaluate a univariate function in B-form.

[mx,nx] = size(x); lx = mx*nx; xs = reshape(x,1,lx);
%  Take apart spline:
[t,a,n,k,d] = spbk(sp);
if lx==0, v = zeros(d,0); return, end

%  Otherwise, augment the knot sequence so that first and last knot each
%  has multiplicity  >= K . (AUGKNT would not be suitable for this since
%  any change in T must be accompanied by a corresponding change in A.)

index = find(diff(t)>0); addl = k-index(1); addr = index(end)-n;
if ( addl>0 || addr>0 )
   npk = n+k; t = t([ones(1,addl) 1:npk npk(ones(1,addr))]);
   a = [zeros(d,addl) a zeros(d,addr)];
   n = n+addl+addr;
end

% For each data point, compute its knot interval:
   [~, index] = histc(-xs,[-inf,-fliplr(t(k+1:n)),inf]);
   NaNx = find(index==0); index = max(n+1-index,k);
if ~isempty(NaNx), index(NaNx) = k; end

% Now, all is ready for the evaluation.
if  k>1  % carry out in lockstep the first spline evaluation algorithm
         % (this requires the following initialization):
   dindex = reshape(repmat(index,d,1),d*lx,1);
   tx =reshape(t(repmat(2-k:k-1,d*lx,1)+repmat(dindex,1,2*(k-1))),d*lx,2*(k-1));
   tx = tx - repmat(reshape(repmat(xs,d,1),d*lx,1),1,2*(k-1));
   dindex = reshape(repmat(d*index,d,1)+repmat((1-d:0).',1,lx),d*lx,1);
   b = repmat(d*(1-k):d:0,d*lx,1)+repmat(dindex,1,k);
   a = a(:); b(:) = a(b);

   % (the following loop is taken from SPRPP)

   for r = 1:k-1
      for i = 1:k-r
         b(:,i) = (tx(:,i+k-1).*b(:,i)-tx(:,i+r-1).*b(:,i+1)) ./ ...
                  (tx(:,i+k-1)    -    tx(:,i+r-1));
      end
   end
   v = reshape(b(:,1),d,lx);
else     % the spline is piecewise constant, hence ...
   v = a(:,index);
   if ~isempty(NaNx), v(:,NaNx) = NaN; end
end

% Finally, zero out all values for points outside the basic interval:
index = find(x<t(1)|x>t(n+k));
if ~isempty(index)
   v(:,index) = zeros(d,length(index));
end
v = reshape(v,d*mx,nx);
end
  
   
function spline = spmk(knots,coefs,sizec)
%based on SPMAK that puts together a spline in B-form.

if nargin<3
 sizec = size(coefs);
end
m = 1;
if iscell(knots)
    m = length(knots);
end

if length(sizec)==m  % coefficients of a scalar-valued function
    sizec = [1 sizec];
end

% convert ND-valued coefficients into vector-valued ones, retaining the
% original size in SIZEVAL, to be temporaily stored in SP.DIM .
sizeval = sizec(1:end-m);
sizec = [prod(sizeval), sizec(end-m+(1:m))];
coefs = reshape(coefs, sizec);

if iscell(knots) % we are putting together a tensor-product spline
    [knots,coefs,k,sizec] = chkt(knots,coefs,sizec);
else            % we are putting together a univariate spline
    [knots,coefs,k,sizec] = chkt({knots},coefs,sizec);
    knots = knots{1};
end

spline.form = 'B-';
spline.knots = knots;
spline.coefs = coefs;
spline.number = sizec(2:end);
spline.order = k;
spline.dim = sizeval;
end

function [knots,coefs,k,sizec] = chkt(knots,coefs,sizec)
%based on CHCKKNT: checks knots and omits trivial B-splines
for j=1:length(sizec)-1
    n = sizec(j+1);
    k(j) = length(knots{j})-n;

       % make sure knot sequence is a row matrix:
    knots{j} = reshape(knots{j},1,n+k(j));
    % throw out trivial B-splines:
    index = find(knots{j}(k(j)+(1:n))-knots{j}(1:n)>0);
    if length(index)<n
        oldn = n;
        n = length(index);
        knots{j} = reshape(knots{j}([index oldn+(1:k(j))]),1,n+k(j));
        coefs = ...
            reshape(coefs, [prod(sizec(1:j)),sizec(j+1),prod(sizec(j+2:end))]);
        sizec(j+1) = n;
        coefs = reshape(coefs(:,index,:),sizec);
    end
end
end


function fprime = fndr(f,dorder)
%based on FNDER:  differentiates a function by differencing the coefficients (see de Boor).
sizeval = f.dim;
if length(sizeval)>1, f.dim = prod(sizeval); end
if nargin<2, dorder=1; end

   [knots,coefs,n,~,d]=spbk(f);
   if iscell(knots)       % the function is multivariate
      m = length(knots);
      sizec = [d,n];% size(coefs);
      for i=m:-1:1
         dsp = fndb(spmk(knots{i},...
            reshape(coefs,prod(sizec(1:m)),sizec(m+1))),dorder(i));
         knots{i} = dsp.knots; sizec(m+1) = dsp.number;
         coefs = reshape(dsp.coefs,sizec); 
         if m>1
            coefs = permute(coefs,[1,m+1,2:m]);
            sizec(2:m+1) = sizec([m+1,2:m]);
         end
      end
      fprime = spmk(knots,coefs,sizec);
   else
      fprime = fndb(f,dorder);
   end
if length(sizeval)>1, fprime.dim = sizeval; end
end


function fprime = fndb(f,dorder)
%FNDERB Differentiate a univariate function in B-form.
[t,a,n,k,d]=spbrk(f);
if k<=dorder
   fprime=spmak(t,zeros(d,n));
elseif dorder<0    % we are to integrate
   error('integration not implemented')
else
   knew=k-dorder;
   for j=k-1:-1:knew
       % here it is: difference the knots and coefficients
      tt=t(j+1+[0:n])-t(1:n+1); z=find(tt>0); nn=length(z);     
      temp=(diff([zeros(1,d);a.'; zeros(1,d)])).';
      a=temp(:,z)./repmat(tt(z)/j,d,1);
      t=[t(z) t(n+2:n+j+1)]; n=nn;
   end
   fprime=spmak(t,a);
end
end


function varargout = spbk(sp)
    varargout = {sp.knots,sp.coefs, sp.number, sp.order, sp.dim};
end


