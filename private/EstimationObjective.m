function [ LogLikelihood, PersistentState, LogObservationLikelihoods ] = EstimationObjective( InputParameters, Options, PersistentState, Smoothing )

    DynamicNu = Options.DynamicNu;
    NoSkewLikelihood = Options.NoSkewLikelihood;
    NoTLikelihood = Options.NoTLikelihood;
    Prior = Options.Prior;
    StationaryDistPeriods = Options.StationaryDistPeriods;
    StationaryDistDrop = Options.StationaryDistDrop;
    StdDevThreshold = Options.StdDevThreshold;
    
    Data = Options.Data;
    Solve = Options.Solve;
    Simulate = Options.Simulate;
    
    ExoCovariance = EstimationOptions.ExoCovariance;
    
    [ T, N ] = size( Data );
    if nargout > 2
        LogObservationLikelihoods = NaN( T, 1 );
    end

    if NoTLikelihood
        Parameters = InputParameters( 1 : ( end - N ) );
        diagLambda = exp( InputParameters( ( end - N + 1 ) : end ) );
        nuoo = Inf;
    elseif DynamicNu
        Parameters = InputParameters( 1 : ( end - N ) );
        diagLambda = exp( InputParameters( ( end - N + 1 ) : end ) );
        nuoo = [];
    else
        Parameters = InputParameters( 1 : ( end - N - 1 ) );
        diagLambda = exp( InputParameters( ( end - N ) : ( end - 1 ) ) );
        nuoo = exp( InputParameters( end ) );
    end
    
    [ PersistentState, StateSteadyState ] = Solve( Parameters, PersistentState );
    
    RootExoCovariance = ObtainEstimateRootCovariance( ExoCovariance, StdDevThreshold );

    OldRNGState = rng( 'default' );
    ShockSequence = RootExoCovariance * randn( size( RootExoCovariance, 2 ), StationaryDistPeriods + StationaryDistDrop );
    rng( OldRNGState );
    
	[ PersistentState, StatDistPoints ] = Simulate( Parameters, PersistentState, ShockSequence, 0 );
    
    StatDistPoints = StatDistPoints( :, ( StationaryDistDrop + 1 ):end );

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

        if isempty( nuoo )
            tauoo_nuoo = lsqnonlin( @( in ) CalibrateMomentsEST( in( 1 ), in( 2 ), MeanStatDist, StateSteadyState, cholVarianceStatDist, sZ3, sZ4 ), [ 2; 10 ], [ -Inf; 4 + eps( 4 ) ], [], optimoptions( @lsqnonlin, 'display', 'off', 'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf ) );
            tauoo = tauoo_nuoo( 1 );
            nuoo = tauoo_nuoo( 2 );
        else
            tauoo = lsqnonlin( @( in ) CalibrateMomentsEST( in( 1 ), nuoo, MeanStatDist, StateSteadyState, cholVarianceStatDist, sZ3, [] ), 2, [], [], optimoptions( @lsqnonlin, 'display', 'off', 'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf ) );
        end
    else
        tauoo = Inf;
        
        if isempty( nuoo )
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

    PriorFunc = str2func( Prior );
    PriorValue = PriorFunc( InputParameters );
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
                nuno = [];
            end
            [ LogObservationLikelihood, xnn, Ssnn, deltasnn, taunn, nunn, wnn, Pnn, deltann, xno, Psno, deltasno, tauno, nuno ] = ...
                KalmanStep( Data( t, : ), xoo, Ssoo, deltasoo, tauoo, nuoo, RootExoCovariance, diagLambda, nuno, PersistentState );
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
            [ LogObservationLikelihood, xnn, Ssnn, deltasnn, taunn, nunn ] = ...
                KalmanStep( Data( t, : ), xoo, Ssoo, deltasoo, tauoo, nuoo, RootExoCovariance, diagLambda, nuno, PersistentState );
            if isempty( xnn )
                error( 'dynareOBC:EstimationEmptyKalmanReturn', 'KalmanStep returned an empty xnn.' );
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
