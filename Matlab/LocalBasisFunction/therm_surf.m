function therm_surf(PT,Y,mask,txt,limits,color)
% function to plot a quantity versus pressure and temperature
% Usage:  therm_surf(PT,Y,mask,txt,limits,color)
% where PT has cells associated with a grid in P and T
%       Y is a matrix (P in columns T in rows) of the surface to plot
%       mask is the same size as Y and has ones where Y should be plotted
%             and NaN where Y values should be supressed in the plot
%       txt is the text to go on the vertical axis
%       limits is a vector defining the plot limits

%colormap gray
GPaflg=1;
P=PT{1};
if(GPaflg)
    P=P/1e3;
end
T=PT{2};

if(nargin<3),
    mask=ones(length(P),length(T));
    txt='  ';
elseif(nargin==3)
    txt='  ';
elseif(nargin==4 & isempty(mask))
    mask=ones(length(P),length(T));
elseif(nargin==5 | nargin==6)
    if(isempty(mask))
      mask=ones(length(P),length(T));
    end
    if(isempty(txt))
        txt=' ';
    end
end
mask0=ones(size(Y));
id=find(Y==0);
mask0(id)=nan;

if nargin<6
    %if(not(isempty(limits)))
    if(exist('limits','var'))
        idl=find(Y<limits(5));
        idh=find(Y>limits(6));
        c=Y;
        c(idl)=limits(5);
        c(idh)=limits(6);
        surf(P,T,mask0'.*mask'.*Y',c')
    else
        surf(P,T,mask0'.*mask'.*Y')
    end
 
 shading 'flat'
else
surf(P,T,mask0'.*mask'.*Y','SpecularExponent',1,...
    'FaceColor',color,...
     'EdgeColor',color,'AlphaData',.3)
end
if(GPaflg)
    xlabel('Pressure (GPa)')
else
    xlabel('Pressure (MPa)')
end
ylabel('Temperature (K)')
zlabel(txt);
if(nargin==5 | nargin==6 & not(isempty(limits)))
    limits(1:2)=limits(1:2)/1e3;
    axis(limits)
end
    
    