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
    Carefully review the options currently supported for these functions.  These are documented in detail
    below.
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

Function options currently supported are listed here.  Unless otherwise stated, are all Boolean values that 
default to _False_.
- _reqM_: if True, concentration (M) is required to calculate the TDV, and the function must include _gPTM_ as one of 
    its input parameters.  See _statevars.evalApparentSpecificHeat_ for an example.
- _reqMWv_: if True, the molecular weight of the solvent is required to calculate the TDV, and the function must include 
    a parameter for the molecular weight of the solvent.  The default name for this parameter is _MWv_, but a different
    parameter name for this value can be specified using  _parmMWv_ when defining the TDV spec in 
    _getSupportedThermodynamicVariables.  See _statevars.evalSolventChemicalPotential_ for an example.
- _reqMWu_: if True, the molecular weight of the solute is required to calculate the TDV, and the function must include  
    a parameter for the molecular weight of the solute.  The default name for this parameter is _MWu_, but a different
    parameter name for this value can be specified using  _parmMWu_ when defining the TDV spec in 
    _getSupportedThermodynamicVariables. See _statevars.evalSoluteChemicalPotential_ for an example.
- _reqGrid_: if True, the PT[M] values must be gridded to properly calculate the TDV, and the function must include a 
    parameter for the gridded PT[M] values.  The default name for this parameter is _gPTM_, but a different parameter
    name for this value can be specified using _parmgrid_ when defining the TDV spec in
    _getSupportedThermodynamicVariables.  See _statevars.evalInternalEnergy_ for an example.
- _reqF_: if True, the TDV requires a conversion factor that gives the volume of solvent in a unit volume of solution
    (which is calculated by evalGibbs._getVolSolventInVolSlnConversion). The function must include a parameter for
    this value, the default name of which is _f_.  A different parameter name can be designated using _parmf_ when 
    defining the TDV spec in _getSupportedThermodynamicVariables.  See _statevars.evalPartialMolarHeatCapacity_ for an 
    example.
- _reqDerivs_: This parameter is not Boolean, but rather a list of derivatives needed to calculate the TDV.  See list of
    supported derivatives in _statevars._getSupportedDerivatives_ (derivatives may be added if necessary). The default 
    value is an empty list, in which case no associated parameter is required for the TDV's calculation function. If any
    derivative is specified, the associated parameter name is _derivs_, unless a different parameter name is specified 
    by the _parmderivs_ parameter when defining the TDV spec in _getSupportedThermodynamicVariables.  See 
    _statevars.evalIsochoricSpecificHeat_ for an example.
- _reqTDV_: This parameter is not Boolean, but rather a list of other TDVs (defined elsewhere in 
    _getSupportedThermodynamicVariables) needed to calculate this TDV.  The default 
    value is an empty list, in which case no associated parameter is required for the TDV's calculation function. If any
    TDV is specified as a dependency, the associated parameter name is _tdv_, unless a different parameter name is 
    specified by the _parmtdv_ parameter when defining the TDV spec in _getSupportedThermodynamicVariables.  See 
    _statevars.evalHelmholtzEnergy_ for an example.  Note that a TDV cannot be dependent on itself.
- _reqSpline_: If True, the spline definition includes additional information required to calculate the TDV.  The 
    eval function must include a parameter to pass in the spline definition.  By default, this parameter should be named
    _gibbsSp_, but a different parameter name can be specified using the _parmspline_ parameter when defining the TDV
    spec in _getSupportedThermodynamicVariables.  See _statevars.evalGibbsEnergy_ for an example.
- _reqPTM_: If True, the pressure/temperature[/concentration] is needed to calculate the TDV, and the eval function must
    include a parameter for passing in the ndarray containing those values.  The default name of that parameter is _PTM_
    but a different parameter name can be specified using the _parmptm_ parameter when defining the TDV spec in 
    _getSupportedThermodynamicVariables_.  See _statevarsevalGibbsEnergy_ for an example.
- _req0M: if True, a value for 0 molal concentration (pure solvent) is required; this is typically used when calculating 
    apparent values. There is no parameter in the eval function associated with this.  Rather, the program will 
    automatically calculate required values for the pure solvent and use them where required.  See 
    _statevars.evalApparentVolume_ for an example.