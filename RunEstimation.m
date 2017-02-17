function [ Parameters, PersistentState ] = RunEstimation( Parameters, Options, PersistentState )

    CorePath = [ fileparts( which( 'RunEstimation' ) ) '/Core/' ];
    addpath( CorePath );

    NumParameters = size( Parameters, 1 );
    NumObservables = size( Options.Data, 1 );
    
    Options = SetDefaultOptions( Options );
    
    InputParameters = [ Parameters; bsxfun( @plus, log( Options.InitialMEStd ), zeros( NumObservables, 1 ) ) ];
        
    EstimatedNu = ~Options.NoTLikelihood && ~Options.DynamicNu;
    if EstimatedNu
        InputParameters( end + 1 ) = log( Options.InitialNu );
    end

    ObjectiveFunction = @( p, s ) EstimationObjective( p, Options, s, false );
    
    [ LogLikelihood, PersistentState ] = ObjectiveFunction( InputParameters, PersistentState );
    PersistentState.InitialRun = false;
    disp( 'Initial log-likelihood:' );
    disp( LogLikelihood );

    LB = Options.LB;
    UB = Options.UB;
    if isempty( LB )
        LB = -Inf( size( Parameters ) );
    end
    if isempty( UB )
        UB = Inf( size( Parameters ) );
    end
    
    OptiLB = [ LB; -Inf( NumObservables + EstimatedNu, 1 ) ];
    OptiUB = [ UB; Inf( NumObservables + EstimatedNu, 1 ) ];
    
    MaximisationFunctions = strsplit( Options.MaximisationFunctions, { ',', ';', '#' } );
    
    for i = 1 : length( MaximisationFunctions )
        FMaxEstimateFunctor = str2func( MaximisationFunctions{ i } );
        [ InputParameters, LogLikelihood, PersistentState ] = FMaxEstimateFunctor( ObjectiveFunction, InputParameters, OptiLB, OptiUB, PersistentState );
    end
    
    disp( 'Final log-likelihood:' );
    disp( LogLikelihood );

    [ LogLikelihood, PersistentState ] = ObjectiveFunction( InputParameters, PersistentState );
    disp( 'Paranoid verification of final log-likelihood:' );
    disp( LogLikelihood );

    Options.PersistentState = PersistentState;

    if Options.SkipStandardErrors
        disp( 'Final parameter estimates:' );
        for i = 1 : NumParameters
            if isempty( Options.ParameterNames ) || length( Options.ParameterNames ) < i
                ParamName = '';
            else
                ParamName = Options.ParameterNames{ i };
            end
            fprintf( '%s:\t\t%#.17g\n', ParamName, InputParameters( i ) );
        end
        fprintf( '\n' );
        disp( 'Final measurement error standard deviation estimates:' );
        for i = 1 : NumObservables
            if isempty( Options.VariableNames ) || length( Options.VariableNames ) < i
                VariableName = '';
            else
                VariableName = Options.VariableNames{ i };
            end
            fprintf( '%s:\t\t%#.17g\n', VariableName, exp( InputParameters( NumParameters + i ) ) );
        end
        if EstimatedNu
            fprintf( '\n' );
            disp( 'Final measurement degrees of freedom parameter:' );
            fprintf( '%s:\t\t%#.17g\n', 'nu', 2 + exp( InputParameters( end ) ) );
        end
    else
        disp( 'Calculating standard errors.' );
        fprintf( '\n' );
        ObservationCount = size( Options.EstimationData, 1 );
        OneOverRootObservationCount = 1 / sqrt( ObservationCount );

        JacobianScoreVector = GetJacobian( @( p ) GetScoreVector( p, Options, PersistentState ), InputParameters, ObservationCount );
        [ ~, TriaJacobianScoreVector ] = qr( JacobianScoreVector * OneOverRootObservationCount, 0 );

        HessianLogLikelihood = GetJacobian( @( p1 ) GetJacobian( @( p2 ) ObjectiveFunction( p2, PersistentState ), p1, 1 )', InputParameters, length( InputParameters ) );
        HessianLogLikelihood = ( 0.5 / ObservationCount ) * ( HessianLogLikelihood + HessianLogLikelihood' );

        RootEstimatedParameterCovarianceMatrix = OneOverRootObservationCount * ( HessianLogLikelihood \ ( TriaJacobianScoreVector' ) );
        EstimatedParameterCovarianceMatrix = RootEstimatedParameterCovarianceMatrix * RootEstimatedParameterCovarianceMatrix';
        Options.EstimatedParameterCovarianceMatrix = EstimatedParameterCovarianceMatrix;
        EstimatedParameterStandardErrors = sqrt( diag( EstimatedParameterCovarianceMatrix ) );

        disp( 'Final parameter estimates:' );
        for i = 1 : NumParameters
            if isempty( Options.ParameterNames ) || length( Options.ParameterNames ) < i
                ParamName = '';
            else
                ParamName = Options.ParameterNames{ i };
            end
            fprintf( '%s:\t\t%#.17g\t\t(%#.17g)\n', ParamName, InputParameters( i ), EstimatedParameterStandardErrors( i ) );
        end
        fprintf( '\n' );
        disp( 'Final measurement error standard deviation estimates:' );
        for i = 1 : NumObservables
            if isempty( Options.VariableNames ) || length( Options.VariableNames ) < i
                VariableName = '';
            else
                VariableName = Options.VariableNames{ i };
            end
            TmpEstimatedParameter = exp( InputParameters( NumParameters + i ) );
            fprintf( '%s:\t\t%#.17g\t\t(%#.17g)\n', VariableName, TmpEstimatedParameter, TmpEstimatedParameter * EstimatedParameterStandardErrors( NumParameters + i ) ); % delta method
        end
        if EstimatedNu
            fprintf( '\n' );
            disp( 'Final measurement degrees of freedom parameter:' );
            TmpEstimatedParameter = exp( InputParameters( end ) );
            fprintf( '%s:\t\t%#.17g\t\t(%#.17g)\n', 'nu', 2 + TmpEstimatedParameter, TmpEstimatedParameter * EstimatedParameterStandardErrors( end ) ); % delta method
        end
    end
    
    Parameters = InputParameters( 1:NumParameters );
    
    rmpath( CorePath );

end
