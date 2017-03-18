function [ CubatureWeights, CubaturePoints, NCubaturePoints ] = GetESTCubaturePoints( xi, Omega, delta, tau, nu, FilterCubatureDegree )

    IntDim = length( xi );
    
    log_tcdf_tau_nu = StudentTLogCDF( tau, nu );
    
    tauTtau = tau * tau;
    OPtauTtauDnu = 1 + tauTtau / nu;
    
    if isfinite( nu )
        nuOnuM1 = nu / ( nu - 1 );
        nuOnuM2 = nu / ( nu - 2 );
    else
        nuOnuM1 = 1;
        nuOnuM2 = 1;
    end
    
    tau2 = tau / realsqrt( nuOnuM2 );
    
    tpdfRatio = exp( StudentTLogPDF( tau, nu ) - log_tcdf_tau_nu );
    
    MedT = -StudentTInvLogCDF( log_tcdf_tau_nu - 0.693147180559945309, nu ); % log( 0.5 ) = -0.693147180559945309
    
    ET1 = nuOnuM1 * OPtauTtauDnu * tpdfRatio;
    ET2 = nuOnuM2 * exp( StudentTLogCDF( tau2, nu - 2 ) - log_tcdf_tau_nu ) - tau * ET1;
    
    lambda = xi + delta * MedT;
    
    if isfinite( nu )
        SamplingCovScale = ( nu + ET2 ) / ( nu - 1 );
    else
        SamplingCovScale = 1;
    end
    
    ET1MMedT = ET1 - MedT;
    
    SamplingCov = SamplingCovScale * Omega + ( max( 0, ET2 - ET1 * ET1 ) + ET1MMedT * ET1MMedT ) * ( delta * delta' );
    
    [ ~, cholSamplingCov ] = NearestSPD( SamplingCov );
    
    [ CubatureWeights, CubaturePoints, NCubaturePoints ] = GetGaussianCubaturePoints( IntDim, FilterCubatureDegree );
    
    CubaturePoints = bsxfun( @plus, lambda, cholSamplingCov' * CubaturePoints );
    
    CubatureWeights = CubatureWeights .* exp( ESTLogPDF( CubaturePoints, xi, Omega, delta, tau, nu ) - MVTNormalLogPDF( CubaturePoints, lambda, cholSamplingCov ) );
    CubatureWeights = CubatureWeights ./ sum( CubatureWeights );
    
end

function [ CubatureWeights, CubaturePoints, NCubaturePoints ] = GetGaussianCubaturePoints( IntDim, FilterCubatureDegree )
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

function log_y = MVTNormalLogPDF( x, mu, cholSigma )
    D = numel( mu );
    log_y = - sum( reallog( abs( diag( cholSigma ) ) ) ) - 0.91893853320467274 * D; % 0.5 * log( 2 * pi ) = 0.91893853320467274
    TIcholSigma_mInnovation = cholSigma' \ bsxfun( @minus, x, mu );
    log_y = log_y - 0.5 * sum( TIcholSigma_mInnovation .* TIcholSigma_mInnovation, 1 );
end
