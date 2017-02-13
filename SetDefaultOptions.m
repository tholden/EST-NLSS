function Options = SetDefaultOptions( Options )
    Options = SetDefaultOption( Options, 'Data', [] );
    Options = SetDefaultOption( Options, 'DynamicNu', false );
    Options = SetDefaultOption( Options, 'ExoCovariance', [] );
    Options = SetDefaultOption( Options, 'FilterCubatureDegree', 0 );
    Options = SetDefaultOption( Options, 'MaximisationFunctions', 'CMAESWrapper,FMinConWrapper' );   
    Options = SetDefaultOption( Options, 'NoSkewLikelihood', false );
    Options = SetDefaultOption( Options, 'NoTLikelihood', false );
    Options = SetDefaultOption( Options, 'Prior', 'FlatPrior' );
    Options = SetDefaultOption( Options, 'Simulate', @( Parameters, PersistentState, ShockSequence )  [] );
    Options = SetDefaultOption( Options, 'Solve', @( Parameters, PersistentState ) [] );
    Options = SetDefaultOption( Options, 'SkipStandardErrors', false );
    Options = SetDefaultOption( Options, 'StationaryDistPeriods', 1000 );
    Options = SetDefaultOption( Options, 'StationaryDistDrop', 100 );
    Options = SetDefaultOption( Options, 'StdDevThreshold', eps ^ 0.375 );    
    Options = orderfields( Options );
end
