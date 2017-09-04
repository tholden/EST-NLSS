function [ log_y, CholOmegaCheck, TIcholOmegaCheck_mInnovation, TIcholOmegaCheck_delta, scaleOmegaNew, scaledeltaNew, tauNew, nuNew ] = ESTLogPDF( x, xi, Omega, delta, tau, nu, OmegaIsCholesky )

    ESTNLSSassert( numel( nu ) == 1, 'ESTNLSS:ESTLogPDF:NuSize', 'ESTLogPDF only supports univariate nu.' );
    ESTNLSSassert( real( nu ) > 0, 'ESTNLSS:ESTLogPDF:NuSign', 'ESTLogPDF requires nu to be strictly positive.' );
    ESTNLSSassert( all( ~isnan( x(:) ) ), 'ESTNLSS:ESTLogPDF:NaNInputX', 'ESTLogPDF was passed a NaN input x.' );
    
    if nargin < 7
        OmegaIsCholesky = false;
    end

    if OmegaIsCholesky
        [ CholOmegaCheck, cholError ] = RealCholeskyUpdate( Omega, delta, '+' );
        if cholError
            Omega = Omega.' * Omega;
            [ ~, CholOmegaCheck ] = NearestSPD( Omega + delta * delta.' );
        end
    else
        [ ~, CholOmegaCheck ] = NearestSPD( Omega + delta * delta.' );
    end

    if coder.target( 'MATLAB' )
        WarningState = warning( 'off', 'MATLAB:nearlySingularMatrix' );
    end
    TIcholOmegaCheck_mInnovation = CholOmegaCheck.' \ bsxfun( @minus, x, xi );
    TIcholOmegaCheck_delta = CholOmegaCheck.' \ delta;
    if coder.target( 'MATLAB' )
        warning( WarningState );
    end

    nx = numel( xi );
    
    if isfinite( real( nu ) )
        scaleOmegaNew = ( nu + sum( TIcholOmegaCheck_mInnovation .* TIcholOmegaCheck_mInnovation, 1 ) ) / ( nu + nx );
    else
        scaleOmegaNew = ones( 1, size( x, 2 ) );
    end
    
    scaledeltaNew = 1 / max( 0, 1 - TIcholOmegaCheck_delta.' * TIcholOmegaCheck_delta );
    tauNew = sqrt( scaledeltaNew ./ scaleOmegaNew ) .* ( sum( bsxfun( @times, TIcholOmegaCheck_delta, TIcholOmegaCheck_mInnovation ), 1 ) + tau );
    nuNew = nu + nx;

    MVTStudentTLogPDF_TIcholOmegaCheck_mInnovation_nu = MVTStudentTLogPDF( TIcholOmegaCheck_mInnovation, nu );
    log_tcdf_tau_nu = ApproxStudentTLogCDF( tau, nu );
    log_tcdf_tauNew_nuNew = ApproxStudentTLogCDF( tauNew, nuNew );
    tcdfDifference = log_tcdf_tauNew_nuNew - log_tcdf_tau_nu;

    log_y = - sum( log( realabs( diag( CholOmegaCheck ) ) ) ) + MVTStudentTLogPDF_TIcholOmegaCheck_mInnovation_nu;

    FiniteCDFDifferenceSelect = isfinite( real( log_tcdf_tau_nu ) ) | isfinite( real( log_tcdf_tauNew_nuNew ) );
    
    log_y( FiniteCDFDifferenceSelect ) = log_y( FiniteCDFDifferenceSelect ) + tcdfDifference( FiniteCDFDifferenceSelect );
        
    ESTNLSSassert( all( ~isnan( log_y(:) ) ), 'ESTNLSS:ESTLogPDF:NaNOutputLogY', 'ESTLogPDF returned a NaN output log_y.' );
    
end

function x = realabs( x )
    SelectNegative = real( x ) < 0;
    x( SelectNegative ) = -x( SelectNegative );
end
