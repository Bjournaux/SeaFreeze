function SF_WPD
%Clinton & Journaux, 2020 
    %This function generates the SeaFreeze Water Phase Diagram. It is calculated
    %by minimizing the differences of Gibbs energy and taking a zero contour. To call the diagram, type:
    %SF_WPD.
    %The data is saved to reduce computing time. See SF_PhaseLines.m for data
    %calculations. 
%Data stored in SeaFreeze WPD.mat
load('WPD','-mat', 'Solid_Solid', 'Melt_Line') 

%Plot Data 
plot(Melt_Line.TP_Ih_water1(2,:),Melt_Line.TP_Ih_water1(1,:),'k-')
hold on
plot(Solid_Solid.TP_Ih_II(2,1:209), Solid_Solid.TP_Ih_II(1,1:209), 'k--');
plot(Solid_Solid.TP_Ih_II(2,209:656), Solid_Solid.TP_Ih_II(1,209:656),'k-');
plot(Solid_Solid.TP_Ih_III(2,:), Solid_Solid.TP_Ih_III(1,:),'k-');
plot(Solid_Solid.TP_II_III(2,:), Solid_Solid.TP_II_III(1,:),'k-');
plot(Solid_Solid.TP_II_V(2,:), Solid_Solid.TP_II_V(1,:),'k-');
plot(Solid_Solid.TP_II_VI(2,506:end), Solid_Solid.TP_II_VI(1,506:end),'k-');
plot(Solid_Solid.TP_II_VI(2,1:506), Solid_Solid.TP_II_VI(1,1:506),'k--');
plot(Solid_Solid.TP_III_V(2,1:19), Solid_Solid.TP_III_V(1,:),'k-');
plot(Melt_Line.TP_III_water1(2,:), Melt_Line.TP_III_water1(1,:),'k-');
plot(Melt_Line.TP_V_water1(2,:), Melt_Line.TP_V_water1(1,:),'k-');
plot(Melt_Line.TP_VI_water1(2,:), Melt_Line.TP_VI_water1(1,:),'k-');
plot(Solid_Solid.TP_V_VI(2,:), Solid_Solid.TP_V_VI(1,:),'k-');

%Add Text
text(50,175,'Ih','FontSize',12)
text(300,300, 'Liquid Water','FontSize',12)
text(400,200,'II','FontSize',12) 
text(260,250, 'III','FontSize',12)
text(525,250, 'V','FontSize',12)
text(1000,250,'VI','FontSize',12)
xlabel('Pressure (MPa)')
ylabel('Temperature (K)') 
%legend('Phase transitions','Metastable extensions','FontSize',12)

Tbound=[0 400];
Pbound=[0 2300];

ylim(Tbound)
xlim(Pbound)