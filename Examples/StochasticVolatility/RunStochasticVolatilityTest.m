Parameters = zeros( 5, 1 );

Parameters( 1 ) = 0;  % log( mu )
Parameters( 2 ) =  3; % log( ( 1 + phi ) / ( 1 - phi ) )
Parameters( 3 ) = -2; % log( omega )
Parameters( 4 ) =  1; % log( ( 1 + rho ) / ( 1 - rho ) )

T = 200;
Drop = 100;
PersistentState = [];

rng( 'default' );

[ PersistentState, TrueEndoSimulation ] = StochasticVolatilitySimulation( Parameters, PersistentState, [], randn( 2, T + Drop ), 0 );

TrueEndoSimulation = TrueEndoSimulation( :, ( Drop + 1 ) : end );
SampleVar = mean( TrueEndoSimulation( 2, : ).^ 2, 2 );
TrueEndoSimulation( 1, : ) = TrueEndoSimulation( 1, : ) - 0.5 * log( SampleVar );
TrueEndoSimulation( 2, : ) = TrueEndoSimulation( 2, : ) ./ SampleVar;
Parameters( 1 ) = Parameters( 1 ) - 0.5 * log( SampleVar );
Data = log( 0.052329478611145268641 + TrueEndoSimulation( 2, : ).^ 2 );

EstimationOptions = struct;

EstimationOptions.AllowTailEvaluations = false;
EstimationOptions.CompileLikelihood = false;
EstimationOptions.Debug = true;
EstimationOptions.DebugMex = true;
EstimationOptions.DynamicNu = true;
EstimationOptions.FilterCubatureDegree = -10;
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
