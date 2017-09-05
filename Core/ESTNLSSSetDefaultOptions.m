function Options = ESTNLSSSetDefaultOptions( Options, Smoothing )
    Options = SetDefaultOption( Options, 'AllowTailEvaluations', false );
    Options = SetDefaultOption( Options, 'CompileLikelihood', false );
    Options = SetDefaultOption( Options, 'Data', [] );
    Options = SetDefaultOption( Options, 'Debug', false );
    Options = SetDefaultOption( Options, 'DebugMex', false );
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
    Options = SetDefaultOption( Options, 'Simulate', @( Parameters, PersistentState, InitialStates, ShockSequence, t ) [] );
    Options = SetDefaultOption( Options, 'Solve', @( Parameters, PersistentState ) [] );
    Options = SetDefaultOption( Options, 'SkipStandardErrors', false );
    Options = SetDefaultOption( Options, 'StationaryDistAccuracy', 10 );
    Options = SetDefaultOption( Options, 'StationaryDistDrop', 0 );
    Options = SetDefaultOption( Options, 'StdDevThreshold', eps ^ 0.375 );    
    Options = SetDefaultOption( Options, 'UB', [] );
    Options = SetDefaultOption( Options, 'MeasurementVariableNames', {} );
    
    Options = SetDefaultOption( Options, 'RootExoCovariance', ObtainEstimateRootCovariance( Options.ExoCovariance, Options.StdDevThreshold ) );
    QMCPoints = HigherOrderSobol( size( Options.RootExoCovariance, 2 ), Options.StationaryDistAccuracy );
    StationaryDistSimulationLength = size( QMCPoints, 2 );
    OldRNGState = rng( 'default' );
    QMCPointsIndices = randperm( StationaryDistSimulationLength );
    rng( OldRNGState );
    Options = SetDefaultOption( Options, 'StationaryDistDraws', QMCPoints( :, QMCPointsIndices ) ); % This procedure is justified by https://arxiv.org/pdf/0807.4858.pdf
    
    assert( Options.StationaryDistDrop < size( Options.StationaryDistDraws, 2 ), 'ESTNLSS:TooHighStationaryDistDrop', 'Options.StationaryDistDrop needs to be lower than the number of points in Options.StationaryDistDraws.' );
    
    if Options.DebugMex
        Options.Debug = true;
    end
    if Options.Debug
        Options.CompileLikelihood = false;
    end
    
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
