function defs = sf_material_defs()
%SF_MATERIAL_DEFS  Centralised material and property metadata for SeaFreeze.
%
%   defs = sf_material_defs()
%
%   Returns a struct with:
%     defs.solid_phases     — ice phase names
%     defs.liquid_phases    — liquid water EOS names
%     defs.nacl_materials   — NaCl(aq) material codes
%     defs.known_materials  — union of the above three
%     defs.base_props       — thermodynamic property names (all materials)
%     defs.shear_props      — shear-modulus-derived property names (solids)
%     defs.mixing_props     — solution mixing property names (NaClaq)
%     defs.MW_H2O           — molar mass of water (kg/mol)

persistent D
if ~isempty(D), defs = D; return; end

D.solid_phases   = {'Ih','II','III','V','VI','VII_X_French'};
D.liquid_phases  = {'water1','water2','water_IAPWS95'};
D.nacl_materials = {'NaClaq','NaClaq_LP','NaClaq_HP','NaClaq_5GPa_2024'};
D.known_materials = [D.solid_phases, D.liquid_phases, D.nacl_materials];

D.base_props   = {'G','S','U','H','A','rho','Cp','Cv','Kt','Kp','Ks', ...
                   'alpha','vel','Js','gamma_Gruneisen','P','T'};
D.shear_props  = {'shear','Vp','Vs'};
D.mixing_props = {'mus','muw','f','m','xs','xw','Va','Cpa','Vm','Vw', ...
                   'Cpm','phi','Vex','aw'};

D.MW_H2O = 0.018015268;  % kg/mol

defs = D;
end
