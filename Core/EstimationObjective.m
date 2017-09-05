function [ LogLikelihood, PersistentState, LogObservationLikelihoods ] = EstimationObjective( EstimatedParameters, Options, PersistentState, Smoothing )

    assert( all( isfinite( EstimatedParameters(:) ) ), 'ESTNLSS:EstimationObjective:NonFiniteInputParameters', 'Non-finite parameters were passed to EstimationObjective.' );
    
    Prior = Options.Prior;
    Solve = Options.Solve;
    Simulate = Options.Simulate;
    
    DynamicNu = Options.DynamicNu;
    SkewLikelihood = ~Options.NoSkewLikelihood;
    NoTLikelihood = Options.NoTLikelihood;
    
    Data = Options.Data;
    
    ExoCovariance = Options.ExoCovariance;
    RootExoCovariance = Options.RootExoCovariance;
    StationaryDistDraws = Options.StationaryDistDraws;
    
    [ N, T ] = size( Data );
    
    if nargout > 2
        LogObservationLikelihoods = NaN( T, 1 );
    end

    if NoTLikelihood
        Parameters = EstimatedParameters( 1 : ( end - N ) );
        diagRootLambda = exp( EstimatedParameters( ( end - N + 1 ) : end ) );
        nuoo = Inf;
    elseif DynamicNu
        Parameters = EstimatedParameters( 1 : ( end - N ) );
        diagRootLambda = exp( EstimatedParameters( ( end - N + 1 ) : end ) );
        nuoo = 0;
    else
        Parameters = EstimatedParameters( 1 : ( end - N - 1 ) );
        diagRootLambda = exp( EstimatedParameters( ( end - N ) : ( end - 1 ) ) );
        nuoo = 5 + exp( EstimatedParameters( end ) );
    end
    
    [ PersistentState.External, StateSteadyState, StateVariableIndices ] = Solve( Parameters, PersistentState.External );
    
    assert( all( isfinite( StateSteadyState(:) ) ), 'ESTNLSS:EstimationObjective:NonFiniteInitialStateSteadyState', 'The solve functor returned a non-finite initial steady-state.' );
    
    ShockSequence = RootExoCovariance * StationaryDistDraws;
    
	[ PersistentState.External, EndoSimulation ] = Simulate( Parameters, PersistentState.External, [], ShockSequence, 0 );
    
    StatDistPoints = EndoSimulation( StateVariableIndices, : );

    StatDistPoints = StatDistPoints( :, ( Options.StationaryDistDrop + 1 ):end );

    assert( all( isfinite( StatDistPoints(:) ) ), 'ESTNLSS:EstimationObjective:NonFiniteStationaryDistSimultation', 'Non-finite values were encountered during the simulation of the stationary distribution.' );
    
    muStatDistPoints = mean( StatDistPoints, 2 );
    SigmaStatDistPoints = cov( StatDistPoints.' );
    [ ~, CholSigmaStatDistPoints ] = NearestSPD( SigmaStatDistPoints );

    if Options.Debug && ~Options.DebugMex
        PersistentState.Internal0 = EstimationObjectiveInternal( StatDistPoints, PersistentState.Internal0, StateSteadyState, DynamicNu, SkewLikelihood, nuoo, muStatDistPoints, CholSigmaStatDistPoints );
    else
        PersistentState.Internal0 = EstimationObjectiveInternal_mex( StatDistPoints, PersistentState.Internal0, StateSteadyState, DynamicNu, SkewLikelihood, nuoo, muStatDistPoints, CholSigmaStatDistPoints );
    end
    
    [ xoo, CholPsoo, deltasoo, tauoo, nuoo ] = GetESTParametersFromVector( PersistentState.Internal0, size( StatDistPoints, 1 ), DynamicNu, SkewLikelihood, nuoo, muStatDistPoints, CholSigmaStatDistPoints );
    
    nuno = nuoo;

    Psoo = CholPsoo * CholPsoo.';

    PriorValue = Prior( EstimatedParameters );
    ScaledPriorValue = PriorValue / T;
    
    if Smoothing
        wnn_ = cell( T, 1 );
        Pnn_ = cell( T, 1 );
        deltann_ = cell( T, 1 );
        taunn_ = cell( T, 1 );
        nunn_ = cell( T, 1 );
        xno_ = cell( T, 1 );
        Psno_ = cell( T, 1 );
        deltasno_ = cell( T, 1 );
        tauno_ = cell( T, 1 );
        nuno_ = cell( T, 1 );
    end
    
    LogLikelihood = 0;   
    
    for t = 1:T
        if DynamicNu
            nuno = 0;
        end
        if Smoothing
            [ PersistentState, LogObservationLikelihood, xnn, Psnn, deltasnn, taunn, nunn, wnn, Pnn, deltann, xno, Psno, deltasno, tauno, nuno ] = ...
                KalmanStep( Data( :, t ), xoo, Psoo, deltasoo, tauoo, nuoo, ExoCovariance, diagRootLambda, nuno, Parameters, Options, PersistentState, StateVariableIndices, t );
            wnn_{ t } = wnn;
            Pnn_{ t } = Pnn;
            deltann_{ t } = deltann;
            taunn_{ t } = taunn;
            nunn_{ t } = nunn;
            xno_{ t } = xno;
            Psno_{ t } = Psno;
            deltasno_{ t } = deltasno;
            tauno_{ t } = tauno;
            nuno_{ t } = nuno;
        else
            [ PersistentState, LogObservationLikelihood, xnn, Psnn, deltasnn, taunn, nunn ] = ...
                KalmanStep( Data( :, t ), xoo, Psoo, deltasoo, tauoo, nuoo, ExoCovariance, diagRootLambda, nuno, Parameters, Options, PersistentState, StateVariableIndices, t );
            if isempty( xnn )
                error( 'ESTNLSS:EmptyKalmanReturn', 'KalmanStep returned an empty xnn.' );
            end
        end
        
        LogObservationLikelihood = LogObservationLikelihood + ScaledPriorValue;
        
        if nargout > 1
            LogObservationLikelihoods( t ) = LogObservationLikelihood;
        end
        
        xoo = xnn;
        Psoo = Psnn;
        deltasoo = deltasnn;
        tauoo = taunn;
        nuoo = nunn;
        
        LogLikelihood = LogLikelihood + LogObservationLikelihood;
    end
    
    if Smoothing
%         SmoothedWs = cell( T, 1 );
%         RootSmoothedWVariances = cell( T, 1 );
%         SmoothedWs{ T } = W;
%         RootSmoothedWVariances{ T } = RootWVariance;
%         
%         for t = ( T - 1 ):-1:1
%             W = FilteredWs{ t } + SmootherGain * ( xnn - PredictedX );
%             SmoothedWs{ t } = W;
%             VarianceTerm1 = RootFilteredWVariances{ t };
%             VarianceTerm2 = SmootherGain * RootPredictedXVariance;
%             VarianceTerm3 = SmootherGain * Ssnn;
%             WVariance = VarianceTerm1 * VarianceTerm1' - VarianceTerm2 * VarianceTerm2' + VarianceTerm3 * VarianceTerm3';
%             RootWVariance = ObtainEstimateRootCovariance( WVariance, 0 );
%             RootSmoothedWVariances{ t } = RootWVariance;
%             Ssnn = RootWVariance( SelectAugStateVariables, : );
%             SmootherGain = SmootherGains{ t };
%             PredictedX = PredictedXs{ t };
%             RootPredictedXVariance = RootPredictedXVariances{ t };
%         end
% 
%         dynareOBC.FilteredWs = FilteredWs;
%         dynareOBC.RootFilteredWVariances = RootFilteredWVariances;
%         dynareOBC.SmoothedWs = SmoothedWs;
%         RootSmoothedWVariances.RootSmoothedWVariances = RootSmoothedWVariances;

    end
        
end
