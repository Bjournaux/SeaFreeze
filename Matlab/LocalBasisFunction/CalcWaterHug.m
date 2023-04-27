function [rhoh,Th,Eh]=CalcWaterHug(sp,ph)
% this function calculates a Hugoniot for water with a variable EOS.
% Usage:  [rhoh,Th,Eh]=CalcWaterHug(sp,ph)
%   where sp is the EOS description and ph is a vector of pressure for the
%   calculation.  rhoh Th and Eh are calculated and returned
% units are MPa, J/kg K and kg/m^3  
% JMB 2015

nh=length(ph);
rhoh=zeros(nh,1);
Th=rhoh;
Eh=rhoh;
[rhoo,~,~,~,~,~,Eo]=fnGval(sp,{0.1,300}); % get initial density and internal energy at STP
Vo=1/rhoo;
Ttrial=300:10:7000;
% it is necessary to determine a valid range of temperatures at each pressure 
mask=mk_mask4Gspline({ph,Ttrial});
%[~,~,~,~,~,~,~,~,mask]=fnIAPWSval(,sp);

% now step through the list of pressures to find the T that satisfies the
% Rankine-Hugoniot energy vs P and V relationship
for i=1:nh
    id=find(mask(i,:)==1);  % valid range for search in T at the specified pressure
    Pin=ph(i);
    Th(i)=fzero(@(x) HugCalc(x),[Ttrial(id(1)) Ttrial(id(end))]); 
    [rhoh(i),~,~,~,~,~,Eh(i)]=fnGval(sp,{ph(i),Th(i)});

end

function f=HugCalc(x)
[rho,~,~,~,~,~,E]=fnGval(sp,{Pin,x});
V=1/rho;
f=(E-Eo)-.5*1e6*Pin*(Vo-V); % pressure in MPa converted to Pa - all other values already in MKS units
end
end
