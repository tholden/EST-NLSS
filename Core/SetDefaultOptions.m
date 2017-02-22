function Options = SetDefaultOptions( Options, Smoothing )
    Options = SetDefaultOption( Options, 'CompileLikelihood', false );
    Options = SetDefaultOption( Options, 'Data', [] );
    Options = SetDefaultOption( Options, 'DynamicNu', false );
    Options = SetDefaultOption( Options, 'ExoCovariance', [] );
    Options = SetDefaultOption( Options, 'FilterCubatureDegree', 0 );
    Options = SetDefaultOption( Options, 'InitialMEStd', 0.0001 );
    Options = SetDefaultOption( Options, 'InitialNu', 20 );
    Options = SetDefaultOption( Options, 'LB', [] );
    Options = SetDefaultOption( Options, 'MaximisationFunctions', { @CMAESWrapper, @FMinConWrapper } );   
    Options = SetDefaultOption( Options, 'NoSkewLikelihood', false );
    Options = SetDefaultOption( Options, 'NoTLikelihood', false );
    Options = SetDefaultOption( Options, 'ParameterNames', {} );
    Options = SetDefaultOption( Options, 'Prior', @FlatPrior );
    Options = SetDefaultOption( Options, 'Simulate', @( Parameters, PersistentState, ShockSequence ) [] );
    Options = SetDefaultOption( Options, 'Solve', @( Parameters, PersistentState ) [] );
    Options = SetDefaultOption( Options, 'SkipStandardErrors', false );
    Options = SetDefaultOption( Options, 'StationaryDistPeriods', 1000 );
    Options = SetDefaultOption( Options, 'StationaryDistDrop', 100 );
    Options = SetDefaultOption( Options, 'StdDevThreshold', eps ^ 0.375 );    
    Options = SetDefaultOption( Options, 'UB', [] );
    Options = SetDefaultOption( Options, 'VariableNames', {} );
    Options = orderfields( Options );
    
    if ~Smoothing && Options.CompileLikelihood
        if ~ischar( Options.Prior )
            Options.Prior = func2str( Options.Prior );
        end
        if ~ischar( Options.Simulate )
            Options.Simulate = func2str( Options.Simulate );
        end
        if ~ischar( Options.Solve )
            Options.Solve = func2str( Options.Solve );
        end
    else
        if ischar( Options.Prior )
            Options.Prior = str2func( Options.Prior );
        end
        if ischar( Options.Simulate )
            Options.Simulate = str2func( Options.Simulate );
        end
        if ischar( Options.Solve )
            Options.Solve = str2func( Options.Solve );
        end
    end
end
