addpath ../..

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

EstimationOptions.ParameterNames = {};
EstimationOptions.VariableNames = {};

EstimationOptions.Data = dynareOBC_.EstimationData;
EstimationOptions.Solve = @EstimationSolution;
EstimationOptions.Simulate = @EstimationSimulation;

NExo = dynareOBC_.OriginalNumVarExo;
EstimationOptions.ExoCovariance = M_.Sigma_e( 1:NExo, 1:NExo );

[ ~, dynareOBC.EstimationParameterSelect ] = ismember( dynareOBC.EstimationParameterNames, cellstr( M_.param_names ) );
LBTemp = dynareOBC.EstimationParameterBounds(1,:)';
UBTemp = dynareOBC.EstimationParameterBounds(2,:)';
LBTemp( ~isfinite( LBTemp ) ) = -Inf;
UBTemp( ~isfinite( UBTemp ) ) = Inf;

EstimationOptions.LB = LBTemp;
EstimationOptions.UB = UBTemp;