function [ resid, xi, delta, CholOmega, Z3, Z4 ] = CalibrateMomentsEST( tau, nu, mu, lambda, CholSigma, sZ3, sZ4 )

    ESTNLSSassert( all( isfinite( mu(:) ) ), 'ESTNLSS:CalibrateMomentsEST:NonFiniteInputMu', 'CalibrateMomentsEST was invoked with a non-finite input mu.' );
    ESTNLSSassert( all( isfinite( lambda(:) ) ), 'ESTNLSS:CalibrateMomentsEST:NonFiniteInputLambda', 'CalibrateMomentsEST was invoked with a non-finite input lambda.' );
    ESTNLSSassert( all( isfinite( CholSigma(:) ) ), 'ESTNLSS:CalibrateMomentsEST:NonFiniteInputCholSigma', 'CalibrateMomentsEST was invoked with a non-finite input CholSigma.' );
    ESTNLSSassert( isempty( sZ3 ) || isfinite( sZ3 ), 'ESTNLSS:CalibrateMomentsEST:NonFiniteInputThirdMoment', 'CalibrateMomentsEST was invoked with a non-finite input third moment.' );
    ESTNLSSassert( isempty( sZ4 ) || isfinite( sZ4 ), 'ESTNLSS:CalibrateMomentsEST:NonFiniteInputFourthMoment', 'CalibrateMomentsEST was invoked with a non-finite input fourth moment.' );
    
    if isfinite( tau )
        tau = max( -10, min( 10, tau ) );
    end
    if isfinite( nu )
        nu = max( 5, min( 1000, nu ) );
    end
    
    log_tcdf_tau_nu = StudentTLogCDF( tau, nu );
    
    resid = zeros( 0, 1 );
    
    ESTNLSSassert( nu > 2, 'ESTNLSS:CalibrateMomentsEST:nuTooSmall', 'nu was too small.' );
    if ~isempty( sZ3 )
        ESTNLSSassert( nu > 3, 'ESTNLSS:CalibrateMomentsEST:nuTooSmall', 'nu was too small.' );
    end
    if ~isempty( sZ4 )
        ESTNLSSassert( nu > 4, 'ESTNLSS:CalibrateMomentsEST:nuTooSmall', 'nu was too small.' );
    end
    
    if log_tcdf_tau_nu == 0
        
        Z3 = 0;
        if isfinite( nu )
            Z4 = ( 3 * nu - 6 ) / ( nu - 4 );
        else
            Z4 = 3;
        end
        if ~isempty( sZ3 )
            resid = [ resid; sZ3 - Z3 ];
        end
        if ~isempty( sZ4 )
            resid = [ resid; sZ4 - Z4 ];
        end
        
        xi = mu;
        delta = zeros( size( mu ) );
        CholOmega = CholSigma;
        if isfinite( nu )
            CholOmega = CholOmega * realsqrt( ( nu - 2 ) / nu );
        end
        
        ESTNLSSassert( all( isfinite( resid(:) ) ), 'ESTNLSS:CalibrateMomentsEST:NonFiniteOutputResid', 'CalibrateMomentsEST returned a non-finite output resid.' );
        ESTNLSSassert( all( isfinite( xi(:) ) ), 'ESTNLSS:CalibrateMomentsEST:NonFiniteOutputXi', 'CalibrateMomentsEST returned a non-finite output xi.' );
        ESTNLSSassert( all( isfinite( delta(:) ) ), 'ESTNLSS:CalibrateMomentsEST:NonFiniteOutputDelta', 'CalibrateMomentsEST returned a non-finite output delta.' );
        ESTNLSSassert( all( isfinite( CholOmega(:) ) ), 'ESTNLSS:CalibrateMomentsEST:NonFiniteOutputResid', 'CalibrateMomentsEST returned a non-finite output CholOmega.' );
        
        return
        
    end

    tauTtau = tau * tau;
    OPtauTtauDnu = 1 + tauTtau / nu;
    
    if isfinite( nu )
        nuOnuM1 = nu / ( nu - 1 );
        nuOnuM2 = nu / ( nu - 2 );
        nuOnuM3 = nu / ( nu - 3 );
    else
        nuOnuM1 = 1;
        nuOnuM2 = 1;
        nuOnuM3 = 1;
    end
    
    tau2 = tau / realsqrt( nuOnuM2 );
    nuTnu = nu * nu;
    
    tpdfRatio = exp( StudentTLogPDF( tau, nu ) - log_tcdf_tau_nu );
    
    MedT = -StudentTInvLogCDF( log_tcdf_tau_nu - 0.693147180559945309, nu ); % log( 0.5 ) = -0.693147180559945309
    
    ET1 = nuOnuM1 * OPtauTtauDnu * tpdfRatio;
    ET2 = nuOnuM2 * exp( StudentTLogCDF( tau2, nu - 2 ) - log_tcdf_tau_nu ) - tau * ET1;
    ET3 = 2 * nuOnuM1 * nuOnuM3 * OPtauTtauDnu * OPtauTtauDnu * tpdfRatio + tauTtau * ET1;
    
    delta = ( mu - lambda ) / ( ET1 - MedT );
    
    xi = mu - delta * ET1;
    delta_deltaT = delta * delta.';
    ET12 = ET1 * ET1;
    
    if isfinite( nu )
        OmegaScaleRatio = ( nu - 1 ) / ( nu + ET2 );
    else
        OmegaScaleRatio = 1;
    end
    
    if ET2 < ET12
        CholOmega = realsqrt( OmegaScaleRatio ) * CholeskyUpdate( CholSigma, realsqrt( ET12 - ET2 ) * delta );
    else
        [ CholOmega, p ] = CholeskyUpdate( CholSigma, realsqrt( ET2 - ET12 ) * delta, '-' );
        if p == 0
            CholOmega = realsqrt( OmegaScaleRatio ) * CholOmega;
        else
            [ ~, CholOmega ] = NearestSPD( OmegaScaleRatio * ( CholSigma.' * CholSigma - ( ET2 - ET12 ) * delta_deltaT ) );
        end
    end
    
    if isempty( sZ3 ) && isempty( sZ4 )
        
        ESTNLSSassert( all( isfinite( resid(:) ) ), 'ESTNLSS:CalibrateMomentsEST:NonFiniteOutputResid', 'CalibrateMomentsEST returned a non-finite output resid.' );
        ESTNLSSassert( all( isfinite( xi(:) ) ), 'ESTNLSS:CalibrateMomentsEST:NonFiniteOutputXi', 'CalibrateMomentsEST returned a non-finite output xi.' );
        ESTNLSSassert( all( isfinite( delta(:) ) ), 'ESTNLSS:CalibrateMomentsEST:NonFiniteOutputDelta', 'CalibrateMomentsEST returned a non-finite output delta.' );
        ESTNLSSassert( all( isfinite( CholOmega(:) ) ), 'ESTNLSS:CalibrateMomentsEST:NonFiniteOutputResid', 'CalibrateMomentsEST returned a non-finite output CholOmega.' );
        
        return
        
    end
    
    cholOmegaCheck = CholeskyUpdate( CholOmega, delta );
    
    deltaT_delta = delta.' * delta;
    deltaT_delta2 = deltaT_delta * deltaT_delta;
    cholOmegaHat_delta = cholOmegaCheck * delta;
    deltaT_OmegaHat_delta = cholOmegaHat_delta' * cholOmegaHat_delta;
    cholSigma_delta = CholSigma * delta;
    deltaT_Sigma_delta = cholSigma_delta' * cholSigma_delta;
    sqrt_deltaT_OmegaHat_delta = realsqrt( deltaT_OmegaHat_delta );
    OmegaHatSigmaRatio = deltaT_OmegaHat_delta / deltaT_Sigma_delta;
    sqrt_delta2OmegaHatRatio = deltaT_delta / sqrt_deltaT_OmegaHat_delta;
    delta2OmegaHatRatio = sqrt_delta2OmegaHatRatio * sqrt_delta2OmegaHatRatio;
    
    omega1 = sqrt_delta2OmegaHatRatio * ET1;
    omega2 = ( deltaT_Sigma_delta + deltaT_delta2 * ET12 ) / deltaT_OmegaHat_delta;
    OMdelta2OmegaHatRatio = 1 - delta2OmegaHatRatio;
    omega3 = 3 * nuOnuM1 * sqrt_delta2OmegaHatRatio * OMdelta2OmegaHatRatio * ( ET1 + ET3 / nu ) + delta2OmegaHatRatio * sqrt_delta2OmegaHatRatio * ET3;

    if ~isempty( sZ3 )
        Z3 = OmegaHatSigmaRatio * realsqrt( OmegaHatSigmaRatio ) * ( omega3 - 3 * omega2 * sqrt_delta2OmegaHatRatio * ET1 + 3 * omega1 * delta2OmegaHatRatio * ET12 - sqrt_delta2OmegaHatRatio * delta2OmegaHatRatio * ET12 * ET1 );
        resid = [ resid; sZ3 - Z3 ];
    end
    
    if ~isempty( sZ4 )
        if isfinite( nu )
            nuOnuM4 = nu / ( nu - 4 );
        else
            nuOnuM4 = 1;
        end
        tau4 = tau / realsqrt( nuOnuM4 );
        if log_tcdf_tau_nu < 0
            ET4 = 3 * nuOnuM2 * nuOnuM4 * exp( StudentTLogCDF( tau4, nu - 4 ) - log_tcdf_tau_nu ) - 1.5 * tau * ET3 + 0.5 * tauTtau * tau * ET1;
        else
            ET4 = 3 * nuOnuM2 * nuOnuM4;
        end
        omega4 = 3 * nuOnuM1 * nuOnuM3 * OMdelta2OmegaHatRatio * OMdelta2OmegaHatRatio * ( 1 + 2 / nu * ET2 + ET4 / nuTnu ) + 6 * nuOnuM1 * delta2OmegaHatRatio * OMdelta2OmegaHatRatio * ( ET2 + ET4 / nu ) + delta2OmegaHatRatio * delta2OmegaHatRatio * ET4;
        Z4 = OmegaHatSigmaRatio * OmegaHatSigmaRatio * ( omega4 - 4 * omega3 * sqrt_delta2OmegaHatRatio * ET1 + 6 * omega2 * delta2OmegaHatRatio * ET12 - 4 * omega1 * sqrt_delta2OmegaHatRatio * delta2OmegaHatRatio * ET12 * ET1 + delta2OmegaHatRatio * delta2OmegaHatRatio * ET12 * ET12 );
        resid = [ resid; sZ4 - Z4 ];
    end
    
    ESTNLSSassert( all( isfinite( resid(:) ) ), 'ESTNLSS:CalibrateMomentsEST:NonFiniteOutputResid', 'CalibrateMomentsEST returned a non-finite output resid.' );
    ESTNLSSassert( all( isfinite( xi(:) ) ), 'ESTNLSS:CalibrateMomentsEST:NonFiniteOutputXi', 'CalibrateMomentsEST returned a non-finite output xi.' );
    ESTNLSSassert( all( isfinite( delta(:) ) ), 'ESTNLSS:CalibrateMomentsEST:NonFiniteOutputDelta', 'CalibrateMomentsEST returned a non-finite output delta.' );
    ESTNLSSassert( all( isfinite( CholOmega(:) ) ), 'ESTNLSS:CalibrateMomentsEST:NonFiniteOutputResid', 'CalibrateMomentsEST returned a non-finite output CholOmega.' );
    
end
