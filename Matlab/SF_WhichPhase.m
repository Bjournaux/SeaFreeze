function out=SF_WhichPhase(PT)
% Version 0.9.1, 
% Calculate which phase is stable at a given PT array
% Returns a numerical value with
% 0 = liquid
% 1 = ice Ih
% 2 = ice II
% 3 = ice III
% 5 = ice V
% 6 = ice VI
% returns NaN for PT values outside the range of all phases
%
%%% Example 
%
% Indentify phase for every 10 MPa from 0 to 2300 MPa and every 1 K from
% 200 to 355 K :
% out = SF_WhichPhase({0:10:2300,200:355})
%

 
load('SeaFreeze_Gibbs.mat')



out_Ih=fnval(G_iceIh,PT');
out_II=fnval(G_iceII,PT');
out_III=fnval(G_iceIII,PT');
out_V=fnval(G_iceV,PT');
out_VI=fnval(G_iceVI,PT');
out_Bollengier=fnval(G_H2O_2GPa_500K,PT');
% out_Brown=fnGval(G_H2O_100GPa_10000K,PT);
% out_IAPWS=fnGval(G_H2O_IAPWS,PT);
[np,nt]=size(out_Ih);

 

 
  for i=1:np*nt
            all_phaseG = [out_Bollengier(i) out_Ih(i) out_II(i) out_III(i) NaN out_V(i) out_VI(i)];
            all_phaseG(find(all_phaseG == 0)) = NaN;
            if all(isnan(all_phaseG))
                out(i) = NaN;
            else
                [Y,I]=min(all_phaseG);
                out(i)=I-1;
            end
  end
   out=reshape(out,np,nt);  
