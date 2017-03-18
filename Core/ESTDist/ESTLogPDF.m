function [ log_y, cholOmegaCheck, TIcholOmegaCheck_mInnovation, TIcholOmegaCheck_delta, scaleOmegaNew, scaledeltaNew, tauNew, nuNew ] = ESTLogPDF( x, xi, Omega, delta, tau, nu )

    [ ~, cholOmegaCheck ] = NearestSPD( Omega + delta * delta' );

    TIcholOmegaCheck_mInnovation = cholOmegaCheck' \ bsxfun( @minus, x, xi );
    TIcholOmegaCheck_delta = cholOmegaCheck' \ delta;

    nx = numel( xi );
    
    if isfinite( nu )
        scaleOmegaNew = ( nu + sum( TIcholOmegaCheck_mInnovation .* TIcholOmegaCheck_mInnovation, 1 ) ) / ( nu + nx );
    else
        scaleOmegaNew = 1;
    end
    
    scaledeltaNew = 1 / ( 1 - TIcholOmegaCheck_delta' * TIcholOmegaCheck_delta );
    tauNew = realsqrt( scaledeltaNew ./ scaleOmegaNew ) .* ( sum( bsxfun( @times, TIcholOmegaCheck_delta, TIcholOmegaCheck_mInnovation ), 1 ) + tau );
    nuNew = nu + nx;

    MVTStudentTLogPDF_TIcholOmegaCheck_mInnovation_nu = MVTStudentTLogPDF( TIcholOmegaCheck_mInnovation, nu );
    log_tcdf_tau_nu = StudentTLogCDF( tau, nu );
    log_tcdf_tauNew_nuNew = StudentTLogCDF( tauNew, nuNew );
    tcdfDifference = log_tcdf_tauNew_nuNew - log_tcdf_tau_nu;

    log_y = - sum( reallog( abs( diag( cholOmegaCheck ) ) ) ) + MVTStudentTLogPDF_TIcholOmegaCheck_mInnovation_nu;

    FiniteCDFDifferenceSelect = isfinite( log_tcdf_tau_nu ) | isfinite( log_tcdf_tauNew_nuNew );
    
    log_y( FiniteCDFDifferenceSelect ) = log_y( FiniteCDFDifferenceSelect ) + tcdfDifference( FiniteCDFDifferenceSelect );
        
end
