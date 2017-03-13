function [ PersistentState, LogObservationLikelihood, xnn, Ssnn, deltasnn, taunn, nunn, wnn, Pnn, deltann, xno, Psno, deltasno, tauno, nuno ] = ...
    KalmanStep( m, xoo, Ssoo, deltasoo, tauoo, nuoo, RootExoVar, diagLambda, nuno, Parameters, Options, PersistentState, StateVariableIndices, t )

    assert( all( isfinite( m(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteInputM', 'KalmanStep was invoked with a non-finite input m.' );
    assert( all( isfinite( xoo(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteInputXoo', 'KalmanStep was invoked with a non-finite input xoo.' );
    assert( all( isfinite( Ssoo(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteInputSsoo', 'KalmanStep was invoked with a non-finite input Ssoo.' );
    assert( all( isfinite( deltasoo(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteInputDeltasoo', 'KalmanStep was invoked with a non-finite input deltasoo.' );
    assert( all( isfinite( RootExoVar(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteInputRootExoVar', 'KalmanStep was invoked with a non-finite input RootExoVar.' );
    assert( all( isfinite( diagLambda(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteInputDiatLambda', 'KalmanStep was invoked with a non-finite input diagLambda.' );
    assert( all( isfinite( Parameters(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteInputParameters', 'KalmanStep was invoked with a non-finite input Parameters.' );

    Simulate = Options.Simulate;
    
%     LogObservationLikelihood = NaN;
%     xnn = [];
%     Ssnn = [];
%     deltasnn = [];
%     taunn = [];
%     nunn = [];
%     wnn = [];
%     Pnn = [];
%     deltann = [];
%     xno = [];
%     Psno = [];
%     deltasno = [];
%     tauno = [];

    NAugState1 = size( Ssoo, 1 );
    NAugState2 = size( Ssoo, 2 );
    NExo1 = size( RootExoVar, 1 );
    NExo2 = size( RootExoVar, 2 );
    
    IntDim = NAugState2 + NExo2 + 2;
    
    tcdf_tauoo_nuoo = StudentTCDF( tauoo, nuoo );
    
    if tcdf_tauoo_nuoo == 1
        IntDim = IntDim - 1;
        tmp_deltasoo = Ssoo \ deltasoo;
        if all( abs( ( Ssoo * tmp_deltasoo - deltasoo ) / max( eps, norm( deltasoo ) ) ) < realsqrt( eps ) )
            % Ssoo * Ssoo' + deltasoo * deltasoo' = Ssoo * Ssoo' + Ssoo * tmp_deltasoo * tmp_deltasoo' * Ssoo' = Ssoo * ( I' * I + tmp_deltasoo * tmp_deltasoo' ) * Ssoo'
            Ssoo = Ssoo * CholeskyUpdate( eye( NAugState2 ), tmp_deltasoo );
        else
            Ssoo = [ Ssoo, deltasoo ];
        end
        deltasoo = zeros( size( deltasoo ) );
    end
    
    if isfinite( nuoo )
        
        [ CubatureWeights, CubaturePoints, NCubaturePoints ] = GetCubaturePoints( IntDim, Options.FilterCubatureDegree );
        PhiN10 = normcdf( CubaturePoints( end, : ) );
        if tcdf_tauoo_nuoo > 0
            N11Scaler = realsqrt( 0.5 * ( nuoo + 1 ) ./ real( gammaincinv( PhiN10, 0.5 * ( nuoo + 1 ), 'upper' ) ) );
        else
            N11Scaler = realsqrt( 0.5 * nuoo ./ real( gammaincinv( PhiN10, 0.5 * nuoo, 'upper' ) ) );
        end
        
        if all( abs( N11Scaler - 1 ) <= realsqrt( eps ) )
            IntDim = IntDim - 1;
            [ CubatureWeights, CubaturePoints, NCubaturePoints ] = GetCubaturePoints( IntDim, Options.FilterCubatureDegree );
            N11Scaler = ones( 1, NCubaturePoints );
        else
            CubaturePoints( end, : ) = [];
        end

    else
        IntDim = IntDim - 1;
        [ CubatureWeights, CubaturePoints, NCubaturePoints ] = GetCubaturePoints( IntDim, Options.FilterCubatureDegree );
        N11Scaler = ones( 1, NCubaturePoints );
    end
    
    if tcdf_tauoo_nuoo < 1 
        PhiN0 = normcdf( CubaturePoints( end, : ) );
        CubaturePoints( end, : ) = [];
        
        ICDFTmp1 = ( 1 - PhiN0 ) * tcdf_tauoo_nuoo;
        ICDFTmp2 = 1 - ICDFTmp1;
        if ICDFTmp2 <= 0.5
            FInvEST = StudentTInvCDF( ICDFTmp2, nuoo );
        else
            FInvEST = -StudentTInvCDF( ICDFTmp1, nuoo );
        end
        
        N11Scaler = N11Scaler .* realsqrt( ( nuoo + FInvEST .^ 2 ) / ( 1 + nuoo ) );
    else
        FInvEST = zeros( 1, NCubaturePoints );
    end

    StateExoPoints = bsxfun( @plus, [ Ssoo * bsxfun( @times, CubaturePoints( 1:NAugState2,: ), N11Scaler ) + bsxfun( @times, deltasoo, FInvEST ); RootExoVar * CubaturePoints( (NAugState2+1):end,: ) ], [ xoo; zeros( NExo1, 1 ) ] );
    
    StatePoints = StateExoPoints( 1:NAugState1, : );
    ExoPoints = StateExoPoints( (NAugState1+1):(NAugState1+NExo1), : );

    assert( all( isfinite( StatePoints( : ) ) ), 'ESTNLSS:KalmanStep:NonFiniteStatePoints', 'Non-finite values were encountered during calculation of the StatePoints in the Kalman step.' );   
    assert( all( isfinite( ExoPoints( : ) ) ), 'ESTNLSS:KalmanStep:NonFiniteExoPoints', 'Non-finite values were encountered during calculation of the ExoPoints in the Kalman step.' );
    
    Observed = find( isfinite( m ) );
    m = m( Observed );
    nm = length( Observed );

    if nm > 0
        [ PersistentState, EndoSimulation, MeasurementSimulation ] = Simulate( Parameters, PersistentState, StatePoints, ExoPoints, t );
        
        assert( all( isfinite( EndoSimulation( : ) ) ), 'ESTNLSS:KalmanStep:NonFiniteSimultation', 'Non-finite values were encountered during simulation in the Kalman step.' );

        NewMeasurementPoints = MeasurementSimulation( Observed, : );
        
        assert( all( isfinite( NewMeasurementPoints(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteMeasurements', 'Non-finite values were encountered during calculation of observation equations in the Kalman step.' );
        
    else
        [ PersistentState, EndoSimulation ] = Simulate( Parameters, PersistentState, StatePoints, ExoPoints, t );
        
        assert( all( isfinite( EndoSimulation( : ) ) ), 'ESTNLSS:KalmanStep:NonFiniteSimultation', 'Non-finite values were encountered during simulation in the Kalman step.' );

        NewMeasurementPoints = zeros( 0, NCubaturePoints );
    end
    
    StdDevThreshold = Options.StdDevThreshold;
    
    NEndo = size( EndoSimulation, 1 );

    wm = [ EndoSimulation; ExoPoints; zeros( nm, NCubaturePoints ); NewMeasurementPoints ];
    
    nwm = size( wm, 1 );
    
    Median_wm = wm( :, 1 );
    
    Mean_wm = sum( bsxfun( @times, wm, CubatureWeights ), 2 );
    ano = bsxfun( @minus, wm, Mean_wm );
    Weighted_ano = bsxfun( @times, ano, CubatureWeights );
    
    Variance_wm = zeros( nwm, nwm );
    ZetaBlock = ( nwm - 2 * nm + 1 ) : nwm;
    Lambda = diag( diagLambda( Observed ) );
    Variance_wm( ZetaBlock, ZetaBlock ) = [ Lambda, Lambda; Lambda, Lambda ];
    
    Variance_wm = Variance_wm + NearestSPD( ano * Weighted_ano' );
    Variance_wm = 0.5 * ( Variance_wm + Variance_wm' );
    [ ~, cholVariance_wm ] = NearestSPD( Variance_wm );
    
    Mean_wmMMedian_wm = Mean_wm - Median_wm;
    cholVariance_wm_Mean_wmMMedian_wm = cholVariance_wm * Mean_wmMMedian_wm;
    cholVariance_wm_Mean_wmMMedian_wm2 = cholVariance_wm_Mean_wmMMedian_wm' * cholVariance_wm_Mean_wmMMedian_wm;
    
    if cholVariance_wm_Mean_wmMMedian_wm2 > eps && ~Options.NoSkewLikelihood
        Zcheck_wm = ( Mean_wmMMedian_wm' * ano ) / realsqrt( cholVariance_wm_Mean_wmMMedian_wm2 );

        meanZcheck_wm = Zcheck_wm * CubatureWeights';
        Zcheck_wm = Zcheck_wm - meanZcheck_wm;
        meanZcheck_wm2 = Zcheck_wm.^2 * CubatureWeights';
        Zcheck_wm = Zcheck_wm / realsqrt( meanZcheck_wm2 );

        sZ3 = Zcheck_wm.^3 * CubatureWeights';
        sZ4 = max( 3, Zcheck_wm.^4 * CubatureWeights' );

        if nuno == 0
            tauno_nuno = LMFnlsq2( @( in ) CalibrateMomentsEST( in( 1 ), 4 + eps( 4 ) + exp( in( 2 ) ), Mean_wm, Median_wm, cholVariance_wm, sZ3, sZ4 ), [ min( 10, tauoo ); log( min( 100, nuoo ) - 4 ) ] );
            tauno = tauno_nuno( 1 );
            nuno = 4 + eps( 4 ) + exp( tauno_nuno( 2 ) );
        else
            tauno = LMFnlsq2( @( in ) CalibrateMomentsEST( in( 1 ), nuno, Mean_wm, Median_wm, cholVariance_wm, sZ3, [] ), min( 10, tauoo ) );
        end
    else
        tauno = Inf;
        
        if nuno == 0
            Zcheck_wm = ano;

            meanZcheck_wm2 = Zcheck_wm.^2 * CubatureWeights';
            Zcheck_wm = bsxfun( @times, Zcheck_wm, 1 ./ realsqrt( meanZcheck_wm2 ) );

            kurtDir = max( 0, Zcheck_wm.^4 * CubatureWeights' - 3 );

            if kurtDir' * kurtDir < eps
                kurtDir = Zcheck_wm.^4 * CubatureWeights';
            end

            kurtDir = kurtDir / norm( kurtDir );

            Zcheck_wm = kurtDir' * Zcheck_wm;

            meanZcheck_wm = Zcheck_wm * CubatureWeights';
            Zcheck_wm = Zcheck_wm - meanZcheck_wm;
            meanZcheck_wm2 = Zcheck_wm.^2 * CubatureWeights';
            Zcheck_wm = Zcheck_wm / realsqrt( meanZcheck_wm2 );

            sZ4 = max( 3, Zcheck_wm.^4 * CubatureWeights' );
            nuno = 4 + 6 / ( sZ4( 1 ) - 3 );
        end
    end
    
    [ ~, wmno, deltaetano, cholPRRQno ] = CalibrateMomentsEST( tauno, nuno, Mean_wm, Median_wm, cholVariance_wm, [], [] );

    assert( NEndo + NExo1 + nm + nm == nwm );
    
    wBlock = 1 : ( NEndo + NExo1 + nm );
    mBlock = ( nwm - nm + 1 ) : nwm;
    
    wno = wmno( wBlock );
    mno = wmno( mBlock );
    
    deltano = deltaetano( wBlock );
    etano = deltaetano( mBlock );
    
    cholPno = cholPRRQno( wBlock, wBlock );
    Pno = cholPno' * cholPno;
    tmpRno = cholPRRQno( wBlock, mBlock );
    Rno = cholPno' * tmpRno;
    tmpQno = cholPRRQno( mBlock, mBlock );
    Qno = tmpRno' * tmpRno + tmpQno' * tmpQno;
    
    xno = wno( StateVariableIndices );
    deltasno = deltano( StateVariableIndices );
    Psno = Pno( StateVariableIndices, StateVariableIndices );
    
    if nm > 0
        cholPnoCheck = CholeskyUpdate( cholPno, deltano );
        
        RnoCheck = Rno + deltano * etano';
        [ ~, cholQnoCheck ] = NearestSPD( Qno + etano * etano' );

        RCheck_IcholQnoCheck = RnoCheck / cholQnoCheck;
        TIcholQnoCheck_mInnovation = cholQnoCheck' \ ( m - mno );
        TIcholQnoCheck_eta = cholQnoCheck' \ etano;
        
        PTildeno = cholPnoCheck * cholPnoCheck' - RCheck_IcholQnoCheck * RCheck_IcholQnoCheck';
        deltaTildeno = deltano - RCheck_IcholQnoCheck * TIcholQnoCheck_eta;
        
        wnn = RCheck_IcholQnoCheck * TIcholQnoCheck_mInnovation;
        if isfinite( nuno )
            scalePnn = ( nuno + TIcholQnoCheck_mInnovation' * TIcholQnoCheck_mInnovation ) / ( nuno + nm );
        else
            scalePnn = 1;
        end
        scaledeltann = 1 / ( 1 - TIcholQnoCheck_eta' * TIcholQnoCheck_eta );
        Pnn = scalePnn * NearestSPD( PTildeno - ( deltaTildeno * deltaTildeno' ) * scaledeltann );
        deltann = realsqrt( scalePnn * scaledeltann ) * deltaTildeno;
        taunn = realsqrt( scaledeltann / scalePnn ) * ( TIcholQnoCheck_eta' * TIcholQnoCheck_mInnovation + tauno );
        nunn = nuno + nm;
        
        [ ~, logMVTStudentTPDF_TIcholQnoCheck_mInnovation_nuno ] = MVTStudentTPDF( TIcholQnoCheck_mInnovation, nuno );
        [ ~, log_tcdf_tauno_nuno ] = StudentTCDF( tauno, nuno );
        [ ~, log_tcdf_taunn_nunn ] = StudentTCDF( taunn, nunn );
        
        LogObservationLikelihood = -sum( log( abs( diag( cholQnoCheck ) ) ) ) + logMVTStudentTPDF_TIcholQnoCheck_mInnovation_nuno;

        if isfinite( log_tcdf_tauno_nuno ) || isfinite( log_tcdf_taunn_nunn )
            tcdfDifference = log_tcdf_taunn_nunn - log_tcdf_tauno_nuno;
            LogObservationLikelihood = LogObservationLikelihood + tcdfDifference;
        end
    else
        wnn = wno;
        Pnn = Pno;
        deltann = deltano;
        taunn = tauno;
        nunn = nuno;
        
        LogObservationLikelihood = 0;
    end
    
    xnn = wnn( StateVariableIndices );
    deltasnn = deltann( StateVariableIndices );
    Psnn = Pnn( StateVariableIndices, StateVariableIndices );
    
    Ssnn = ObtainEstimateRootCovariance( Psnn, StdDevThreshold );
       
    % [ PersistentState, LogObservationLikelihood, xnn, Ssnn, deltasnn, taunn, nunn, wnn, Pnn, deltann, xno, Psno, deltasno, tauno, nuno ] = ...

    assert( isfinite( LogObservationLikelihood ), 'ESTNLSS:KalmanStep:NonFiniteOutputLogObservationLikelihood', 'KalmanStep returned a non-finite output log observation likelihood.' );
    assert( all( isfinite( xnn(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteOutputXnn', 'KalmanStep returned a non-finite output xnn.' );
    assert( all( isfinite( Ssnn(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteOutputSsnn', 'KalmanStep returned a non-finite output Ssnn.' );
    assert( all( isfinite( deltasnn(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteOutputDeltasnn', 'KalmanStep returned a non-finite output deltasnn.' );
    assert( all( isfinite( wnn(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteOutputWnn', 'KalmanStep returned a non-finite output wnn.' );
    assert( all( isfinite( Pnn(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteOutputPnn', 'KalmanStep returned a non-finite output Pnn.' );
    assert( all( isfinite( deltann(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteOutputDeltann', 'KalmanStep returned a non-finite output deltann.' );
    assert( all( isfinite( xno(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteOutputXno', 'KalmanStep returned a non-finite output xno.' );
    assert( all( isfinite( Psno(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteOutputPsno', 'KalmanStep returned a non-finite output Psno.' );
    assert( all( isfinite( deltasno(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteOutputDeltasno', 'KalmanStep returned a non-finite output deltasno.' );


end

function [ CubatureWeights, CubaturePoints, NCubaturePoints ] = GetCubaturePoints( IntDim, FilterCubatureDegree )
    if FilterCubatureDegree > 0
        CubatureOrder = ceil( 0.5 * ( FilterCubatureDegree - 1 ) );
        [ CubatureWeights, CubaturePoints, NCubaturePoints ] = fwtpts( IntDim, CubatureOrder );
    else
        NCubaturePoints = 2 * IntDim + 1;
        wTemp = 0.5 * realsqrt( 2 * NCubaturePoints );
        CubaturePoints = [ zeros( IntDim, 1 ), wTemp * eye( IntDim ), -wTemp * eye( IntDim ) ];
        CubatureWeights = 1 / NCubaturePoints;
    end
end
