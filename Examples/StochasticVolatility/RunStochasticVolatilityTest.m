Parameters = zeros( 5, 1 );

Parameters( 1 ) = -4; % logit( mu )
Parameters( 2 ) =  3; % log( ( 1 + phi ) / ( 1 - phi ) )
Parameters( 3 ) = -2; % log( omega )
Parameters( 4 ) =  1; % log( ( 1 + rho ) / ( 1 - rho ) )
Parameters( 5 ) =  0; % log( max_sigma )

T = 200;
Drop = 100;
PersistentState = [];

rng( 'default' );

[ PersistentState, TrueEndoSimulation, Data ] = StochasticVolatilitySimulation( Parameters, PersistentState, [], randn( 2, T + Drop ), 0 );

TrueEndoSimulation = TrueEndoSimulation( :, ( Drop + 1 ) : end );
Data = Data( :, ( Drop + 1 ) : end );

EstimationOptions = struct;

EstimationOptions.AllowTailEvaluations = false;
EstimationOptions.CompileLikelihood = false;
EstimationOptions.Debug = false;
EstimationOptions.DebugMex = false;
EstimationOptions.DynamicNu = true;
EstimationOptions.FilterCubatureDegree = 19;
EstimationOptions.MaximisationFunctions = 'FMinConWrapper';
EstimationOptions.NoSkewLikelihood = false;
EstimationOptions.NoTLikelihood = false;
EstimationOptions.Prior = '';
EstimationOptions.SkipStandardErrors = false;
EstimationOptions.StationaryDistAccuracy = 10;
EstimationOptions.StdDevThreshold = eps;

EstimationOptions.ParameterNames = { 'log( mu )', 'log( ( 1 + phi ) / ( 1 - phi ) )', 'log( omega )', 'log( ( 1 + rho ) / ( 1 - rho ) )' };
EstimationOptions.VariableNames = { 'log( sigma )', 'e' };

EstimationOptions.Data = Data;
EstimationOptions.Solve = @StochasticVolatilitySolution;
EstimationOptions.Simulate = @StochasticVolatilitySimulation;

EstimationOptions.ExoCovariance = eye( 2 );

EstimationOptions.LB = -Inf( 4, 1 );
EstimationOptions.UB = Inf( 4, 1 );

addpath ../..

[ EstimatedParameters, PersistentState ] = RunEstimation( Parameters, EstimationOptions, PersistentState );
