function [ LogLikelihood, PersistentState, LogObservationLikelihoods ] = EstimationObjective( EstimatedParameters, Options, PersistentState, Smoothing )

    Prior = Options.Prior;
    Solve = Options.Solve;
    Simulate = Options.Simulate;
    
    DynamicNu = Options.DynamicNu;
    NoSkewLikelihood = Options.NoSkewLikelihood;
    NoTLikelihood = Options.NoTLikelihood;
    StationaryDistPeriods = Options.StationaryDistPeriods;
    StationaryDistDrop = Options.StationaryDistDrop;
    StdDevThreshold = Options.StdDevThreshold;
    
    Data = Options.Data;
    
    ExoCovariance = Options.ExoCovariance;
    
    [ N, T ] = size( Data );
    
    if nargout > 2
        LogObservationLikelihoods = NaN( T, 1 );
    end

    if NoTLikelihood
        Parameters = EstimatedParameters( 1 : ( end - N ) );
        diagLambda = exp( 2 * EstimatedParameters( ( end - N + 1 ) : end ) );
        nuoo = Inf;
    elseif DynamicNu
        Parameters = EstimatedParameters( 1 : ( end - N ) );
        diagLambda = exp( 2 * EstimatedParameters( ( end - N + 1 ) : end ) );
        nuoo = 0;
    else
        Parameters = EstimatedParameters( 1 : ( end - N - 1 ) );
        diagLambda = exp( 2 * EstimatedParameters( ( end - N ) : ( end - 1 ) ) );
        nuoo = exp( EstimatedParameters( end ) );
    end
    
    [ PersistentState, StateSteadyState, StateVariableIndices ] = Solve( Parameters, PersistentState );
    
    RootExoCovariance = ObtainEstimateRootCovariance( ExoCovariance, StdDevThreshold );

    OldRNGState = rng( 'default' );
    ShockSequence = RootExoCovariance * randn( size( RootExoCovariance, 2 ), StationaryDistPeriods + StationaryDistDrop );
    rng( OldRNGState );
    
	[ PersistentState, EndoSimulation ] = Simulate( Parameters, PersistentState, [], ShockSequence, 0 );
    
    StatDistPoints = EndoSimulation( StateVariableIndices, : );

    StatDistPoints = StatDistPoints( :, ( StationaryDistDrop + 1 ):end );

    if any( ~isfinite( StatDistPoints(:) ) )
        error( 'ESTNLSS:NonFiniteStationaryDistSimultation', 'Non-finite values were encountered during the simulation of the stationary distribution.' );
    end

    MeanStatDist = mean( StatDistPoints, 2 );
    DeMeanedStatDistPoints = bsxfun( @minus, StatDistPoints, MeanStatDist );
    [ ~, cholVarianceStatDist ] = NearestSPD( cov( StatDistPoints' ) );

    MeanStatDistMMedianStatDist = MeanStatDist - StateSteadyState;
    cholVarianceStatDist_MeanStatDistMMedianStatDist = cholVarianceStatDist * MeanStatDistMMedianStatDist;
    cholVarianceStatDist_MeanStatDistMMedianStatDist2 = cholVarianceStatDist_MeanStatDistMMedianStatDist' * cholVarianceStatDist_MeanStatDistMMedianStatDist;
    
    nuno = nuoo;
    
    if cholVarianceStatDist_MeanStatDistMMedianStatDist2 > eps && ~NoSkewLikelihood
        ZcheckStatDist = ( MeanStatDistMMedianStatDist' * DeMeanedStatDistPoints ) / sqrt( cholVarianceStatDist_MeanStatDistMMedianStatDist2 );

        sZ3 = skewness( ZcheckStatDist, 0 );
        sZ4 = max( 3, kurtosis( ZcheckStatDist, 0 ) );

        if nuoo == 0
            tauoo_nuoo = LMFnlsq2( @( in ) CalibrateMomentsEST( in( 1 ), 4 + eps( 4 ) + exp( in( 2 ) ), MeanStatDist, StateSteadyState, cholVarianceStatDist, sZ3, sZ4 ), [ 2; 2 ] );
            tauoo = tauoo_nuoo( 1 );
            nuoo = 4 + eps( 4 ) + exp( tauoo_nuoo( 2 ) );
        else
            tauoo = LMFnlsq2( @( in ) CalibrateMomentsEST( in( 1 ), nuoo, MeanStatDist, StateSteadyState, cholVarianceStatDist, sZ3, [] ), 2 );
        end
    else
        tauoo = Inf;
        
        if nuoo == 0
            kurtDir = max( 0, kurtosis( DeMeanedStatDistPoints, 0, 2 ) - 3 );

            if kurtDir' * kurtDir < eps
                kurtDir = kurtosis( DeMeanedStatDistPoints, 0, 2 );
            end

            kurtDir = kurtDir / norm( kurtDir );

            ZcheckStatDist = kurtDir' * DeMeanedStatDistPoints;

            sZ4 = max( 3, kurtosis( ZcheckStatDist, 0 ) );
            nuoo = 4 + 6 / ( sZ4 - 3 );
        end
    end
    
    [ ~, xoo, deltasoo, cholPsoo ] = CalibrateMomentsEST( tauoo, nuoo, MeanStatDist, StateSteadyState, cholVarianceStatDist, [], [] );
    
    Psoo = cholPsoo * cholPsoo';
    Ssoo = ObtainEstimateRootCovariance( Psoo, StdDevThreshold );

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
        if Smoothing
            if DynamicNu
                nuno = 0;
            end
            [ PersistentState, LogObservationLikelihood, xnn, Ssnn, deltasnn, taunn, nunn, wnn, Pnn, deltann, xno, Psno, deltasno, tauno, nuno ] = ...
                KalmanStep( Data( :, t ), xoo, Ssoo, deltasoo, tauoo, nuoo, RootExoCovariance, diagLambda, nuno, Parameters, Options, PersistentState, StateVariableIndices, t );
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
            [ PersistentState, LogObservationLikelihood, xnn, Ssnn, deltasnn, taunn, nunn ] = ...
                KalmanStep( Data( :, t ), xoo, Ssoo, deltasoo, tauoo, nuoo, RootExoCovariance, diagLambda, nuno, Parameters, Options, PersistentState, StateVariableIndices, t );
            if isempty( xnn )
                error( 'ESTNLSS:EmptyKalmanReturn', 'KalmanStep returned an empty xnn.' );
            end
        end
        
        LogObservationLikelihood = LogObservationLikelihood + ScaledPriorValue;
        
        if nargout > 1
            LogObservationLikelihoods( t ) = LogObservationLikelihood;
        end
        
        xoo = xnn;
        Ssoo = Ssnn;
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
