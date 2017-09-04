function [ PersistentState, LogObservationLikelihood, xnn, Psnn, deltasnn, taunn, nunn, wnn, Pnn, deltann, xno, Psno, deltasno, tauno, nuno ] = ...
    KalmanStep( m, xoo, Psoo, deltasoo, tauoo, nuoo, ExoVar, diagRootLambda, nuno, Parameters, Options, PersistentState, StateVariableIndices, t )

    assert( all( isfinite( m(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteInputM', 'KalmanStep was invoked with a non-finite input m.' );
    assert( all( isfinite( xoo(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteInputXoo', 'KalmanStep was invoked with a non-finite input xoo.' );
    assert( all( isfinite( Psoo(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteInputSsoo', 'KalmanStep was invoked with a non-finite input Ssoo.' );
    assert( all( isfinite( deltasoo(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteInputDeltasoo', 'KalmanStep was invoked with a non-finite input deltasoo.' );
    assert( all( isfinite( ExoVar(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteInputRootExoVar', 'KalmanStep was invoked with a non-finite input RootExoVar.' );
    assert( all( isfinite( diagRootLambda(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteInputDiatLambda', 'KalmanStep was invoked with a non-finite input diagRootLambda.' );
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

    NAugState1 = size( Psoo, 1 );
    NAugState2 = size( Psoo, 2 );
    NExo1 = size( ExoVar, 1 );
    NExo2 = size( ExoVar, 2 );
    
    ZerosNExo1 = zeros( NExo1, 1 );
    
    [ CubatureWeights, StateExoPoints, NCubaturePoints ] = GetESTCubaturePoints( [ xoo; ZerosNExo1 ], [ Psoo, zeros( NAugState1, NExo2 ); zeros( NExo1, NAugState2 ), ExoVar ], [ deltasoo; ZerosNExo1 ], tauoo, nuoo, Options.FilterCubatureDegree, Options.StdDevThreshold, Options.AllowTailEvaluations );
    
    StatePoints = StateExoPoints( 1:NAugState1, : );
    ExoPoints = StateExoPoints( (NAugState1+1):(NAugState1+NExo1), : );

    assert( all( isfinite( StatePoints( : ) ) ), 'ESTNLSS:KalmanStep:NonFiniteStatePoints', 'Non-finite values were encountered during calculation of the StatePoints in the Kalman step.' );   
    assert( all( isfinite( ExoPoints( : ) ) ), 'ESTNLSS:KalmanStep:NonFiniteExoPoints', 'Non-finite values were encountered during calculation of the ExoPoints in the Kalman step.' );
    
    Observed = find( isfinite( m ) );
    m = m( Observed );
    nm = length( Observed );

    if nm > 0
        [ PersistentState.External, EndoSimulation, MeasurementSimulation ] = Simulate( Parameters, PersistentState.External, StatePoints, ExoPoints, t );
        
        assert( all( isfinite( EndoSimulation( : ) ) ), 'ESTNLSS:KalmanStep:NonFiniteSimultation', 'Non-finite values were encountered during simulation in the Kalman step.' );

        NewMeasurementPoints = MeasurementSimulation( Observed, : );
        
        assert( all( isfinite( NewMeasurementPoints(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteMeasurements', 'Non-finite values were encountered during calculation of observation equations in the Kalman step.' );
        
    else
        [ PersistentState.External, EndoSimulation ] = Simulate( Parameters, PersistentState.External, StatePoints, ExoPoints, t );
        
        assert( all( isfinite( EndoSimulation( : ) ) ), 'ESTNLSS:KalmanStep:NonFiniteSimultation', 'Non-finite values were encountered during simulation in the Kalman step.' );

        NewMeasurementPoints = zeros( 0, NCubaturePoints );
    end

    % wm = [ EndoSimulation; ExoPoints; zeros( nm, NCubaturePoints ); NewMeasurementPoints ];
    
    wmRed = [ EndoSimulation; ExoPoints; NewMeasurementPoints ];
    
    DynamicNu = nuno == 0;
    SkewLikelihood = ~Options.NoSkewLikelihood;
    
    nwmRed = size( wmRed, 1 );
    
    ParamDim = ( nwmRed + 1 ) * SkewLikelihood + DynamicNu;
    
    PersistentState.Internal( ( end + 1 ):ParamDim, : ) = NaN;
    
    mu_wmRed = sum( bsxfun( @times, wmRed, CubatureWeights ), 2 );
    
    Demeaned_wmRed = bsxfun( @minus, wmRed, mu_wmRed );
    Weighted_Demeaned_wmRed = bsxfun( @times, Demeaned_wmRed, CubatureWeights );

    Sigma_wmRed = Demeaned_wmRed * Weighted_Demeaned_wmRed.';
    [ ~, CholSigma_wmRed ] = NearestSPD( Sigma_wmRed );
    
    if Options.Debug && ~Options.DebugMex
        PersistentState.Internal( 1 : ParamDim, t ) = KalmanStepInternal( wmRed, CubatureWeights,  PersistentState.Internal( 1 : ParamDim, t ), DynamicNu, SkewLikelihood, nuno, mu_wmRed, CholSigma_wmRed );
    else
        PersistentState.Internal( 1 : ParamDim, t ) = KalmanStepInternal_mex( wmRed, CubatureWeights,  PersistentState.Internal( 1 : ParamDim, t ), DynamicNu, SkewLikelihood, nuno, mu_wmRed, CholSigma_wmRed );
    end
    
    [ wmno, CholPRRQno, deltaetano, tauno, nuno ] = GetESTParametersFromVector( PersistentState.Internal( 1 : ParamDim, t ), nwmRed, DynamicNu, SkewLikelihood, nuno, mu_wmRed, CholSigma_wmRed );
    
    NEndo = size( EndoSimulation, 1 );
    
    assert( NEndo + NExo1 + nm == nwmRed );
    
    % add back the nm elements of zeta which are missing here
    
    lFirstBlock = NEndo + NExo1;
    FirstBlock = 1 : lFirstBlock;
    SecondBlock = ( lFirstBlock + 1 ) : ( lFirstBlock + nm );
    
    MECholPRRQnoComponent = sqrt( GetOmegaScaleRatio( tauno, nuno ) ) * diag( diagRootLambda( Observed ) );
    
    wmno = [ wmno( FirstBlock ); zeros( nm, 1 ); wmno( SecondBlock ) ];
    CholPRRQno = [ CholPRRQno( FirstBlock, FirstBlock ), zeros( lFirstBlock, nm ), CholPRRQno( FirstBlock, SecondBlock ); zeros( nm, lFirstBlock ), MECholPRRQnoComponent, MECholPRRQnoComponent; zeros( nm, nwmRed ), CholPRRQno( SecondBlock, SecondBlock ) ];
    deltaetano = [ deltaetano( FirstBlock ); zeros( nm, 1 ); deltaetano( SecondBlock ) ];  

    nwm = lFirstBlock + 2 * nm;
    
    assert( size( wmno, 1 ) == nwm );

    wBlock = 1 : ( lFirstBlock + nm );
    mBlock = ( lFirstBlock + nm + 1 ) : nwm;
    
    wno = wmno( wBlock );
    mno = wmno( mBlock );
    
    deltano = deltaetano( wBlock );
    etano = deltaetano( mBlock );
    
    CholPno = CholPRRQno( wBlock, wBlock );
    Pno = CholPno.' * CholPno;
    tmpRno = CholPRRQno( wBlock, mBlock );
    Rno = CholPno.' * tmpRno;
    tmpQno = CholPRRQno( mBlock, mBlock );
    Qno = tmpRno.' * tmpRno + tmpQno.' * tmpQno;
    
    xno = wno( StateVariableIndices );
    deltasno = deltano( StateVariableIndices );
    Psno = Pno( StateVariableIndices, StateVariableIndices );
    
    if nm > 0
        CholPnoCheck = CholeskyUpdate( CholPno, deltano );
        
        RnoCheck = Rno + deltano * etano.';
        
        [ LogObservationLikelihood, CholQnoCheck, TICholQnoCheck_mInnovation, TICholQnoCheck_eta, scalePnn, scaledeltann, taunn, nunn ] = ESTLogPDF( m, mno, Qno, etano, tauno, nuno );

        RCheck_ICholQnoCheck = RnoCheck / CholQnoCheck;
        
        PTildeno = CholPnoCheck.' * CholPnoCheck - RCheck_ICholQnoCheck * RCheck_ICholQnoCheck.';
        deltaTildeno = deltano - RCheck_ICholQnoCheck * TICholQnoCheck_eta;
        
        wnn = RCheck_ICholQnoCheck * TICholQnoCheck_mInnovation;

        Pnn = scalePnn * NearestSPD( PTildeno - ( deltaTildeno * deltaTildeno.' ) * scaledeltann );
        deltann = realsqrt( scalePnn * scaledeltann ) * deltaTildeno;        
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
    
    % [ PersistentState, LogObservationLikelihood, xnn, Ssnn, deltasnn, taunn, nunn, wnn, Pnn, deltann, xno, Psno, deltasno, tauno, nuno ] = ...

    assert( isfinite( LogObservationLikelihood ), 'ESTNLSS:KalmanStep:NonFiniteOutputLogObservationLikelihood', 'KalmanStep returned a non-finite output log observation likelihood.' );
    assert( all( isfinite( xnn(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteOutputXnn', 'KalmanStep returned a non-finite output xnn.' );
    assert( all( isfinite( Psnn(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteOutputSsnn', 'KalmanStep returned a non-finite output Ssnn.' );
    assert( all( isfinite( deltasnn(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteOutputDeltasnn', 'KalmanStep returned a non-finite output deltasnn.' );
    assert( all( isfinite( wnn(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteOutputWnn', 'KalmanStep returned a non-finite output wnn.' );
    assert( all( isfinite( Pnn(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteOutputPnn', 'KalmanStep returned a non-finite output Pnn.' );
    assert( all( isfinite( deltann(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteOutputDeltann', 'KalmanStep returned a non-finite output deltann.' );
    assert( all( isfinite( xno(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteOutputXno', 'KalmanStep returned a non-finite output xno.' );
    assert( all( isfinite( Psno(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteOutputPsno', 'KalmanStep returned a non-finite output Psno.' );
    assert( all( isfinite( deltasno(:) ) ), 'ESTNLSS:KalmanStep:NonFiniteOutputDeltasno', 'KalmanStep returned a non-finite output deltasno.' );

end

function OmegaScaleRatio = GetOmegaScaleRatio( tau, nu )
    log_tcdf_tauno_nuno = StudentTLogCDF( tau, nu );
    
    if isfinite( nu )
        nuOnuM1 = nu / ( nu - 1 );
        nuOnuM2 = nu / ( nu - 2 );
    else
        nuOnuM1 = 1;
        nuOnuM2 = 1;
    end
    
    tau2 = tau / realsqrt( nuOnuM2 );
    
    tauTtau = tau * tau;
    
    if tauTtau < Inf
        ET1 = nuOnuM1 * ( 1 + tauTtau / nu ) * exp( StudentTLogPDF( tau, nu ) - log_tcdf_tauno_nuno );
        ET2 = nuOnuM2 * exp( StudentTLogCDF( tau2, nu - 2 ) - log_tcdf_tauno_nuno ) - tau * ET1;
    elseif tau >= 0
        ET2 = nuOnuM2;
    else
        ET2 = Inf;
    end
    
    if isfinite( nu )
        OmegaScaleRatio = ( nu - 1 ) / ( nu + ET2 );
    else
        OmegaScaleRatio = 1;
    end
end
