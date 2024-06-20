**Note for developers:**


Logic in this module should only deal with generic MatLab B-spline evaluation.  
There should be no logic related to thermodynamics -- any such logic should instead be added to the 
``lbftd`` module.  The ``loadSpline`` function output will include any additional parameters included 
in the spline that are unused or unrecognized by generic spline evaluation.