Currently supports only pure substances or single-solute solutions

__**Warning: units must be as specified here because some conversions are hardcoded.**__
With the exception of pressure, units are SI. 
 - Pressure is in MPa (rather than the SI standard Pa).
 - Temperature is in K.
 - Concentration (molality) is in mol kg<sup>-1</sup>.
 
Supported thermodynamic variables (TDVs) are of two types: those that require concentration (M) to calculate and those 
that rely only on pressure (P) and temperature (T).

These rely only on P and T:
- _G_:  Gibbs energy in J kg<sup>-1</sup>
- _rho_: density in kg m<sup>-3</sup>
- _vel_: sound speed in m s<sup>-1</sup>
- _Cp_: isobaric specific heat in J kg<sup>-1</sup> K<sup>-1</sup>
- _Cv_: isochoric specific heat in J kg<sup>-1</sup> K<sup>-1</sup>
- _alpha_: thermal expansivity in K<sup>-1</sup>
- _U_: internal energy in J kg<sup>-1</sup>
- _A_: Helmholtz energy in J kg<sup>-1</sup>
- _H_: enthalpy in J kg<sup>-1</sup>
- _S_: entropy in J kg<sup>-1</sup> K<sup>-1</sup>
- _Kt_: isothermal bulk modulus in MPa
- _Kp_: pressure derivatives of isothermal bulk modulus (dimensionless)
- _Ks_: isotropic bulk modulus in MPa
- _V_: unit volume in m<sup>3</sup> kg<sup>-1</sup>

These rely on P,T, and M, and some also require non-zero molecular weights (for solvent and solute) to calculate.
- _mus_: solute chemical potential in J mol<sup>-1</sup>
- _muw_: solvent chemical potential in J mol<sup>-1</sup>
- _Vm_: partial molar volume in m<sup>3</sup> mol<sup>-1</sup>
- _Cpm_: partial molar heat capacity in J kg<sup>-1</sup> K<sup>-1</sup> mol<sup>-1</sup>
- _Cpa_: apparent heat capacity J kg<sup>-1</sup> K<sup>-1</sup> mol<sup>-1</sup>
- _Va_: apparent volume m<sup>3</sup> mol<sup>-1</sup>


See demo folder for example usage.  


__**For developers**__:

To add a new thermodynamic variable (TDV), all of the following should be done.  This list may not be comprehensive.
-  New variables cannot be named PTM, P, T, or M, as those symbols are reserved for the input. 
    Read the comments for statevars.getSupportedMeasures and statevars._getTDVSpec.
- Create a short function in statevars to calculate the measure based on other values
    such as gibbsSp, PTM, gPTM, derivs, tdvout, etc.  
    The procedure should be named with 'eval' + the FULL name for the measure - NOT the symbol / tdv flag.  
    Record the symbol / tdv flag as the return value in comments for the function.
    Add only parameters required to directly calculate the measure.
    Be consistent with parameter names used in other functions or use the parm* parameters of 
    statevars._getTDVSpec.
    If you end up with an as-yet unused parameter, add it to statevars._getTDVSpec (defaulting to OFF)
    AND to the evaluation section of evalGibbs.evalSolutionGibbs.
- Add the measure spec to statevars.getSupportedThermodynamicVariables.
    When the comments say DIRECTLY, they mean only consider something a requirement if it is used in
    the function built in the previous step.
    Dependencies (including nested dependencies) will be handled by the reqDerivs and reqTDV parameters.
    See statevars functions _addTDVDependencies and expandTDVspec, as well as evalGibbs.getDerivatives, for detail.
- Update this README with the name of the measure *and its units*.
    Be sure to add it to the correct section of the comments (PT vs PTM spline, other parameters required, etc)
    or create a new section if one is warranted.
- Add tests to make sure that the TDV spec still expands properly and that the values are calculated correctly.  
    The latter may require recalculating a spline in MatLab and comparing it with the output from that platform. 