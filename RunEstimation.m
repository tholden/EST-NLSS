function [ EstimatedParameters, EstimatedParameterCovarianceMatrix, PersistentState ] = RunEstimation( Parameters, Options, PersistentState )

    CorePath = [ fileparts( which( 'RunEstimation' ) ) '/Core/' ];
    addpath( CorePath );
    
    try

        NumParameters = size( Parameters, 1 );
        NumObservables = size( Options.Data, 1 );

        Options = SetDefaultOptions( Options, false );

        EstimatedParameters = [ Parameters; bsxfun( @plus, log( Options.InitialMEStd ), zeros( NumObservables, 1 ) ) ];

        EstimatedNu = ~Options.NoTLikelihood && ~Options.DynamicNu;
        if EstimatedNu
            EstimatedParameters( end + 1 ) = log( Options.InitialNu );
        end

        if Options.CompileLikelihood
            cfg = coder.config( 'mex' );
            cfg.EnableMemcpy = false;
            cfg.InitFltsAndDblsToZero = false;
            cfg.ConstantInputs = 'Remove';
            cfg.MATLABSourceComments = true;
            cfg.GenerateReport = true;
            cfg.ConstantFoldingTimeout = 2147483647;
            cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
            cfg.SaturateOnIntegerOverflow = false;
            cfg.EnableAutoExtrinsicCalls = false;
            cfg.InlineThreshold = 2147483647;
            cfg.InlineThresholdMax = 2147483647;
            cfg.InlineStackLimit = 2147483647;
            cfg.StackUsageMax = 16777216;
            cfg.IntegrityChecks = false;
            cfg.ResponsivenessChecks = false;
            cfg.ExtrinsicCalls = false;
            cfg.EchoExpressions = false;
            cfg.GlobalDataSyncMethod = 'NoSync';

            TmpOptions = Options;
            TmpOptions = rmfield( TmpOptions, 'ParameterNames' );
            TmpOptions = rmfield( TmpOptions, 'MeasurementVariableNames' );

            ARGS = cell( 4, 1 );
            ARGS{1} = coder.typeof( EstimatedParameters );
            ARGS{2} = coder.Constant( TmpOptions );
            ARGS{3} = coder.typeof( PersistentState );
            ARGS{4} = coder.Constant( false ); %#ok<NASGU>

            EstimationObjectiveText = fileread( 'EstimationObjective.m' );
            EstimationObjectiveLines = strsplit( EstimationObjectiveText, { '\n', '\r' } );
            EstimationObjectiveLines( 2:4 ) = [];
            EstimationObjectiveText = strjoin( EstimationObjectiveLines, '\n' );
            EstimationObjectiveText = strrep( EstimationObjectiveText, 'Prior', Options.Prior );
            EstimationObjectiveText = strrep( EstimationObjectiveText, 'Solve', Options.Solve );
            EstimationObjectiveText = strrep( EstimationObjectiveText, 'Simulate', Options.Simulate );
            EstimationObjectiveText = strrep( EstimationObjectiveText, 'EstimationObjective', 'ESTNLSSTempEstimationObjective' );
            EstimationObjectiveText = strrep( EstimationObjectiveText, 'KalmanStep', 'ESTNLSSTempKalmanStep' );

            TempEstimationObjectiveFID = fopen( 'ESTNLSSTempEstimationObjective.m', 'w' );
            fprintf( TempEstimationObjectiveFID, '%s', EstimationObjectiveText );
            fclose( TempEstimationObjectiveFID );

            KalmanStepText = fileread( 'KalmanStep.m' );
            KalmanStepLines = strsplit( KalmanStepText, { '\n', '\r' } );
            KalmanStepLines( 3 ) = [];
            KalmanStepText = strjoin( KalmanStepLines, '\n' );
            KalmanStepText = strrep( KalmanStepText, 'Simulate', Options.Simulate );
            KalmanStepText = strrep( KalmanStepText, 'KalmanStep', 'ESTNLSSTempKalmanStep' );

            TempKalmanStepFID = fopen( 'ESTNLSSTempKalmanStep.m', 'w' );
            fprintf( TempKalmanStepFID, '%s', KalmanStepText );
            fclose( TempKalmanStepFID );

            codegen -config cfg ESTNLSSTempEstimationObjective -args ARGS -o ESTNLSSTempEstimationObjectiveMex;
            rehash;
            ObjectiveFunction = @ESTNLSSTempEstimationObjectiveMex;
        else
            ObjectiveFunction = @( p, s ) EstimationObjective( p, Options, s, false );
        end

        [ LogLikelihood, PersistentState ] = ObjectiveFunction( EstimatedParameters, PersistentState );
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

        MaximisationFunctions = Options.MaximisationFunctions;
        
        for i = 1 : length( MaximisationFunctions )
            FMaxEstimateFunctor = MaximisationFunctions{ i };
            if ischar( FMaxEstimateFunctor )
                FMaxEstimateFunctor = str2func( FMaxEstimateFunctor );
            end
            [ EstimatedParameters, LogLikelihood, PersistentState ] = FMaxEstimateFunctor( ObjectiveFunction, EstimatedParameters, OptiLB, OptiUB, PersistentState );
        end

        disp( 'Final log-likelihood:' );
        disp( LogLikelihood );

        [ LogLikelihood, PersistentState ] = ObjectiveFunction( EstimatedParameters, PersistentState );
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
                fprintf( '%s:\t\t%#.17g\n', ParamName, EstimatedParameters( i ) );
            end
            fprintf( '\n' );
            disp( 'Final measurement error standard deviation estimates:' );
            for i = 1 : NumObservables
                if isempty( Options.MeasurementVariableNames ) || length( Options.MeasurementVariableNames ) < i
                    VariableName = '';
                else
                    VariableName = Options.MeasurementVariableNames{ i };
                end
                fprintf( '%s:\t\t%#.17g\n', VariableName, exp( EstimatedParameters( NumParameters + i ) ) );
            end
            if EstimatedNu
                fprintf( '\n' );
                disp( 'Final measurement degrees of freedom parameter:' );
                fprintf( '%s:\t\t%#.17g\n', 'nu', 2 + exp( EstimatedParameters( end ) ) );
            end
            EstimatedParameterCovarianceMatrix = [];
        else
            disp( 'Calculating standard errors.' );
            fprintf( '\n' );
            ObservationCount = size( Options.EstimationData, 1 );
            OneOverRootObservationCount = 1 / sqrt( ObservationCount );

            JacobianScoreVector = GetJacobian( @( p ) GetScoreVector( p, PersistentState, ObjectiveFunction ), EstimatedParameters, ObservationCount );
            [ ~, TriaJacobianScoreVector ] = qr( JacobianScoreVector * OneOverRootObservationCount, 0 );

            HessianLogLikelihood = GetJacobian( @( p1 ) GetJacobian( @( p2 ) ObjectiveFunction( p2, PersistentState ), p1, 1 )', EstimatedParameters, length( EstimatedParameters ) );
            HessianLogLikelihood = ( 0.5 / ObservationCount ) * ( HessianLogLikelihood + HessianLogLikelihood' );

            RootEstimatedParameterCovarianceMatrix = OneOverRootObservationCount * ( HessianLogLikelihood \ ( TriaJacobianScoreVector' ) );
            EstimatedParameterCovarianceMatrix = RootEstimatedParameterCovarianceMatrix * RootEstimatedParameterCovarianceMatrix';
            EstimatedParameterStandardErrors = sqrt( diag( EstimatedParameterCovarianceMatrix ) );

            disp( 'Final parameter estimates:' );
            for i = 1 : NumParameters
                if isempty( Options.ParameterNames ) || length( Options.ParameterNames ) < i
                    ParamName = '';
                else
                    ParamName = Options.ParameterNames{ i };
                end
                fprintf( '%s:\t\t%#.17g\t\t(%#.17g)\n', ParamName, EstimatedParameters( i ), EstimatedParameterStandardErrors( i ) );
            end
            fprintf( '\n' );
            disp( 'Final measurement error standard deviation estimates:' );
            for i = 1 : NumObservables
                if isempty( Options.MeasurementVariableNames ) || length( Options.MeasurementVariableNames ) < i
                    VariableName = '';
                else
                    VariableName = Options.MeasurementVariableNames{ i };
                end
                TmpEstimatedParameter = exp( EstimatedParameters( NumParameters + i ) );
                fprintf( '%s:\t\t%#.17g\t\t(%#.17g)\n', VariableName, TmpEstimatedParameter, TmpEstimatedParameter * EstimatedParameterStandardErrors( NumParameters + i ) ); % delta method
            end
            if EstimatedNu
                fprintf( '\n' );
                disp( 'Final measurement degrees of freedom parameter:' );
                TmpEstimatedParameter = exp( EstimatedParameters( end ) );
                fprintf( '%s:\t\t%#.17g\t\t(%#.17g)\n', 'nu', 2 + TmpEstimatedParameter, TmpEstimatedParameter * EstimatedParameterStandardErrors( end ) ); % delta method
            end
        end
        
        Error = [];
    
    catch Error
    end
    
    rmpath( CorePath );
    
    if ~isempty( Error )
        rethrow( Error );
    end

end
