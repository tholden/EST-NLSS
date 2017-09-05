function [ xi, CholOmega, delta, tau, nu ] = GetESTParametersFromVector( p, n, DynamicNu, SkewLikelihood, nu, mu, CholSigma )

    if SkewLikelihood
        delta = p( ( 1 : n ).' );
        tau = p( n + 1 );
    else
        delta = zeros( n, 1 );
        tau = Inf;
    end
    if DynamicNu
        nu = 5 + exp( p( end ) );
    end

    if isfinite( nu )
        nuOnuM1 = nu / ( nu - 1 );
        nuOnuM2 = nu / ( nu - 2 );
    else
        nuOnuM1 = 1;
        nuOnuM2 = 1;
    end
    
    log_tcdf_tau_nu = ApproxStudentTLogCDF( tau, nu );
    
    if log_tcdf_tau_nu == 0
        
        ET1 = 0;
        ET2 = nuOnuM2;
        
    else
        
        tauTtau = tau * tau;
        OPtauTtauDnu = 1 + tauTtau / nu;

        tau2 = tau / sqrt( nuOnuM2 );

        tpdfRatio = exp( StudentTLogPDF( tau, nu ) - log_tcdf_tau_nu );

        ET1 = nuOnuM1 * OPtauTtauDnu * tpdfRatio;
        ET2 = nuOnuM2 * exp( ApproxStudentTLogCDF( tau2, nu - 2 ) - log_tcdf_tau_nu ) - tau * ET1;
        
    end
    
    xi = mu - delta * ET1;
    
    delta_deltaT = delta * delta.';
    ET12 = ET1 * ET1;
    
    if isfinite( nu )
        OmegaScaleRatio = ( nu - 1 ) / ( nu + ET2 );
    else
        OmegaScaleRatio = 1;
    end
    
    if ET2 < ET12
        [ CholOmega, CholError ] = RealCholeskyUpdate( CholSigma, sqrt( ET12 - ET2 ) * delta, '+' );
    else
        [ CholOmega, CholError ] = RealCholeskyUpdate( CholSigma, sqrt( ET2 - ET12 ) * delta, '-' );
    end
    
    if CholError
        [ ~, CholOmega ] = NearestSPD( OmegaScaleRatio * ( CholSigma.' * CholSigma - ( ET2 - ET12 ) * delta_deltaT ) );
    else
        CholOmega = sqrt( OmegaScaleRatio ) * CholOmega;
    end    
    
end
