Parameters = zeros( 4, 1 );

Parameters( 1 ) = -4; % log( mu )
Parameters( 2 ) =  3; % log( ( 1 + phi ) / ( 1 - phi ) )
Parameters( 3 ) = -2; % log( omega )
Parameters( 4 ) =  1; % log( ( 1 + rho ) / ( 1 - rho ) )

T = 1000;
Drop = 100;
PersistentState = [];
[ PersistentState, TrueEndoSimulation, Data ] = StochasticVolatilitySimulation( Parameters, PersistentState, [], randn( 2, T + Drop ), 0 );

TrueEndoSimulation = TrueEndoSimulation( :, ( Drop + 1 ) : end );
Data = Data( :, ( Drop + 1 ) : end );

EstimationOptions = struct;

EstimationOptions.DynamicNu = true;
EstimationOptions.FilterCubatureDegree = 9;
EstimationOptions.MaximisationFunctions = 'FMinConWrapper';
EstimationOptions.NoSkewLikelihood = false;
EstimationOptions.NoTLikelihood = false;
EstimationOptions.Prior = @FlatPrior;
EstimationOptions.SkipStandardErrors = false;
EstimationOptions.StationaryDistPeriods = 1000;
EstimationOptions.StationaryDistDrop = 100;
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
