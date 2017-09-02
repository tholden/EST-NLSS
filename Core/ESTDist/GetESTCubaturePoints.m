function [ CubatureWeights, CubaturePoints, NCubaturePoints, ET1, MedT ] = GetESTCubaturePoints( xi, Omega, delta, tau, nu, FilterCubatureDegree, StdDevThreshold, AllowTailEvaluations )

    log_tcdf_tau_nu = StudentTLogCDF( tau, nu );
    
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
        ET1 = nuOnuM1 * ( 1 + tauTtau / nu ) * exp( StudentTLogPDF( tau, nu ) - log_tcdf_tau_nu );
        ET2 = nuOnuM2 * exp( StudentTLogCDF( tau2, nu - 2 ) - log_tcdf_tau_nu ) - tau * ET1;
        VarT = max( 0, ET2 - ET1 * ET1 );
    elseif tau >= 0
        ET1 = 0;
        ET2 = nuOnuM2;
        VarT = nuOnuM2;
    else
        ET1 = Inf;
        ET2 = Inf;
        VarT = 0;
    end
    
    MedT = -StudentTInvLogCDF( log_tcdf_tau_nu - 0.693147180559945309, nu ); % log( 0.5 ) = -0.693147180559945309
    
    lambda = xi + delta * MedT;
    
    if isfinite( nu )
        SamplingCovScale = ( nu + ET2 ) / ( nu - 1 );
    else
        SamplingCovScale = 1;
    end
    
    ET1MMedT = ET1 - MedT;
    
    SamplingCov = SamplingCovScale * Omega + ( VarT + ET1MMedT * ET1MMedT ) * ( delta * delta' );
    
    RootSamplingCov = ObtainEstimateRootCovariance( SamplingCov, StdDevThreshold );
    
    IntDim = size( RootSamplingCov, 2 );
    
    [ CubatureWeights, UncorrelatedCubaturePoints, NCubaturePoints ] = GetGaussianCubaturePoints( IntDim, FilterCubatureDegree );
    
    if AllowTailEvaluations
        SamplingNu = nu + IntDim - 1; % Lowest nu such that the ratio of desired pdf / sampling pdf is bounded below
    else
        SamplingNu = Inf;
    end
    
    ESTNLSSassert( SamplingNu > 2, 'ESTNLSS:GetESTCubaturePoints:SamplingNuTooSmall', 'GetESTCubaturePoints requires that nu + IntDim - 1 is greater than 2.' );
    
    if isfinite( SamplingNu )
        RootSamplingCov = RootSamplingCov * sqrt( ( SamplingNu - 2 ) / SamplingNu );
        UncorrelatedCubaturePointSigns = sign( UncorrelatedCubaturePoints );
        UncorrelatedCubaturePoints = -abs( UncorrelatedCubaturePoints );
        UncorrelatedCubaturePoints = UncorrelatedCubaturePointSigns .* StudentTInvLogCDF( StudentTLogCDF( UncorrelatedCubaturePoints, Inf ), SamplingNu );
    end
    
    CubaturePoints = bsxfun( @plus, lambda, RootSamplingCov * UncorrelatedCubaturePoints );
    
    CubatureWeights = CubatureWeights .* exp( ESTLogPDF( CubaturePoints, xi, Omega, delta, tau, nu ) - sum( StudentTLogPDF( UncorrelatedCubaturePoints, SamplingNu ) ) ); % no need to include the term in the determinant of RootSamplingCov as it washes out in the next line
    CubatureWeights = CubatureWeights ./ sum( CubatureWeights );
    
end
