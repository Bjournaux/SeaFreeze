function out=SeaFreeze(PT,material)
% Version 0.9.2 ; Journaux et al. 2020
% Calculate thermodynamic quantities for water or ices polymorphs 
% (Ih, II, III, V and VI). Needs the SeaFreeze_Gibbs.mat library containing the
% Gibbs Local Basis Function parametrization to run
% Reference publication : Journaux et al., (2020)
% usage:
%    out=SeaFreeze(PT,'material')
% where:
%     out is a structure containing all output quantities (SI units):
%           - G (J/kg) Gibbs Energy
%           - S (J/K/kg) Entropy
%           - U (J/kg) Internal Energy
%           - H (J/kg) Enthalpy
%           - rho(kg/m^3) Density
%           - Cp (J/kg/K) Specific heat capacity at constant pressure
%           - Cv (J/kg/K) Specific heat capacity at constant volume
%           - Kt (MPa) Isothermal bulk modulus
%           - Kp (dimensionless) Pressure derivative of the Isothermal bulk modulus
%           - Ks (MPa) Isoentropic bulk modulus
%           - alpha (/K)  Thermal expansivity
%           - shear (MPa) Shear modulus
%           - Vp (m/s) P wave velocity
%           - Vs (m/s) S wave velocity
%           - vel (m/s) bulk sound speed
%     PT is a structure (gridded output) or array (scatter output)
%          containing pressure-temperature points (MPa and Kelvin)
%     NaN values returned when out of parametrization boundaries.
%
%    Material defines which ice or water to use.  Possibilities:
%         Ih for ice Ih (Feistel and Wagner, 2006)
%         II for ice II (Journaux et al. 2020)
%         III for ice III (Journaux et al. 2020)
%         V for ice V (Journaux et al. 2020)
%         VI for ice VI (Journaux et al. 2020)
%         water1 for Bollengier et al. (2019) LBF extending to 500 K and 2300 MPa
%         water2 for the modified EOS in Brown 2018 extending to 100 GPa
%         water_IAPWS95 for IAPWS95 water (Wagner and Pruß, 2002)
%
%%%% !! Important remarks : !!
%       The ices Gibbs parametrizations are optimized to be used with
%       'water1' Gibbs LBF from Bollengier et al. (2019), specially for phase
%       equilibrium calculation. Using other water parametrization wil lead
%       to incorect melting curves.
%       'water2' (Brown 2018) and 'water_IAPWS95' (IAPWS95) parametrization
%       are provided for HP extention (up to 100 GPa) and comparison only. 
%       The authors recommend the use of 'water1' (Bollengier et al. 2019)
%       for any application in the 200-355 K range and up to 2300 MPa.
%
%%%% Examples  :
%  - Single point for ice VI at 900 MPa and 255 K
%          out=SeaFreeze([900,255],'VI')
%  - Array for ice V every 2 MPa from 400 to 500 MPa and every 0.5 K from
%  220 to 250 K
%          out=SeaFreeze({400:2:500,240:0.5:250},'VI')
%
%
%%%% References
%
%   Bollengier, Brown and Shaw (2019) J. Chem. Phys. 151; doi: 10.1063/1.5097179
%   Brown (2018) Fluid Phase Equilibria 463, pp. 18-31
%   Feistel and Wagner (2006), J. Phys. Chem. Ref. Data 35, pp. 1021-1047
%   Journaux et al., (2020) JGR: Planets, 125, e2019JE006176.
%   Wagner and Pruß (2002), J. Phys. Chem. Ref. Data 31, pp. 387-535

test_toolbox = license('test','curve_fitting_toolbox');



% spline tools from curve_fitting_toolbox are faster but require the 
% toolbox to be installed. An alternative (OpenGval) based on the  sp_val function
% is provided to allow use of LBF representations without installation of 
% the optional Curevefitting Toolbox. It is as fast for gridded PT inputs
% and around 5 times slower for a list of inputs
if test_toolbox == 1;   
    switch material
        case 'Ih'
            load('SeaFreeze_Gibbs.mat', 'G_iceIh');
            out=fnGval(G_iceIh,PT);
            shear_mod=[3.1 -0.00462 0 -0.00657 1000 273.15];

        case 'II'
            load('SeaFreeze_Gibbs.mat', 'G_iceII');
             out=fnGval(G_iceII,PT);
             shear_mod=[4.1 0.0175 0 -0.014 1100 273];
            
        case 'III'
            load('SeaFreeze_Gibbs.mat', 'G_iceIII');
             out=fnGval(G_iceIII,PT);
             shear_mod=[2.57 0.0175 0 -0.014 1100 273];
        case 'V'
            load('SeaFreeze_Gibbs.mat', 'G_iceV');
             out=fnGval(G_iceV,PT);
             shear_mod=[2.57 0.0175 0 -0.014 1100 273];
        case 'VI'
            load('SeaFreeze_Gibbs.mat', 'G_iceVI');
             out=fnGval(G_iceVI,PT);
            shear_mod=[2.57 0.0175 0 -0.014 1100 273];

        case 'water1'
            load('SeaFreeze_Gibbs.mat', 'G_H2O_2GPa_500K');
             out=fnGval(G_H2O_2GPa_500K,PT);
        case 'water2'
            load('SeaFreeze_Gibbs.mat', 'G_H2O_100GPa_10000K');
             out=fnGval(G_H2O_100GPa_10000K,PT);
        case 'water_IAPWS95'
            load('SeaFreeze_Gibbs.mat', 'G_H2O_IAPWS');
            out=fnGval(G_H2O_IAPWS,PT);
    end
else
    switch material
        case 'Ih'
            load('SeaFreeze_Gibbs.mat', 'G_iceIh');
            out=OpenGval(G_iceIh,PT);
            shear_mod=[3.1 -0.00462 0 -0.00657 1000 273.15];
        
        case 'II'
            load('SeaFreeze_Gibbs.mat', 'G_iceII');
             out=OpenGval(G_iceII,PT);
             shear_mod=[4.1 0.0175 0 -0.014 1100 273];
            
        case 'III'
            load('SeaFreeze_Gibbs.mat', 'G_iceIII');
             out=OpenGval(G_iceIII,PT);
             shear_mod=[2.57 0.0175 0 -0.014 1100 273];
        case 'V'
            load('SeaFreeze_Gibbs.mat', 'G_iceV');
             out=OpenGval(G_iceV,PT);
             shear_mod=[2.57 0.0175 0 -0.014 1100 273];
        case 'VI'
            load('SeaFreeze_Gibbs.mat', 'G_iceVI');
             out=OpenGval(G_iceVI,PT);
            shear_mod=[2.57 0.0175 0 -0.014 1100 273];

        case 'water1'
            load('SeaFreeze_Gibbs.mat', 'G_H2O_2GPa_500K');
             out=OpenGval(G_H2O_2GPa_500K,PT);
        case 'water2'
            load('SeaFreeze_Gibbs.mat', 'G_H2O_100GPa_10000K');
             out=OpenGval(G_H2O_100GPa_10000K,PT);
        case 'water_IAPWS95'
            load('SeaFreeze_Gibbs.mat', 'G_H2O_IAPWS');
            out=OpenGval(G_H2O_IAPWS,PT);
    end
end
    
    
    
if iscell(PT) 
    [~,Tm]=ndgrid(PT{1},PT{2});
else
    Tm=PT(:,2);
end
if(exist('shear_mod'))
    shear=shear_mod(1)+shear_mod(2)*(out.rho-shear_mod(5))+shear_mod(3)*(out.rho-shear_mod(5)).^2+shear_mod(4)*(Tm-shear_mod(6));
    out.Vp=1e3*sqrt((out.Ks/1e3+4/3*shear)./out.rho/1e-3);
    out.Vs=1e3*sqrt(shear./out.rho/1e-3);
    out.shear=shear*1e3; % convert to MPa to be consistent with units
end

