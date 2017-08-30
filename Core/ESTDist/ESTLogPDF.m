function [ log_y, cholOmegaCheck, TIcholOmegaCheck_mInnovation, TIcholOmegaCheck_delta, scaleOmegaNew, scaledeltaNew, tauNew, nuNew ] = ESTLogPDF( x, xi, Omega, delta, tau, nu )

    assert( numel( nu ) == 1, 'ESTNLSS:ESTLogPDF:NuSize', 'ESTLogPDF only supports univariate nu.' );
    assert( real( nu ) > 0, 'ESTNLSS:ESTLogPDF:NuSign', 'ESTLogPDF requires nu to be strictly positive.' );
    assert( all( ~isnan( x(:) ) ), 'ESTNLSS:ESTLogPDF:NaNInputX', 'ESTLogPDF was passed a NaN input x.' );

    [ ~, cholOmegaCheck ] = NearestSPD( Omega + delta * delta' );

    TIcholOmegaCheck_mInnovation = cholOmegaCheck' \ bsxfun( @minus, x, xi );
    TIcholOmegaCheck_delta = cholOmegaCheck' \ delta;

    nx = numel( xi );
    
    if isfinite( real( nu ) )
        scaleOmegaNew = ( nu + sum( TIcholOmegaCheck_mInnovation .* TIcholOmegaCheck_mInnovation, 1 ) ) / ( nu + nx );
    else
        scaleOmegaNew = 1;
    end
    
    scaledeltaNew = 1 / ( 1 - TIcholOmegaCheck_delta' * TIcholOmegaCheck_delta );
    tauNew = sqrt( scaledeltaNew ./ scaleOmegaNew ) .* ( sum( bsxfun( @times, TIcholOmegaCheck_delta, TIcholOmegaCheck_mInnovation ), 1 ) + tau );
    nuNew = nu + nx;

    MVTStudentTLogPDF_TIcholOmegaCheck_mInnovation_nu = MVTStudentTLogPDF( TIcholOmegaCheck_mInnovation, nu );
    log_tcdf_tau_nu = ApproxStudentTLogCDF( tau, nu );
    log_tcdf_tauNew_nuNew = ApproxStudentTLogCDF( tauNew, nuNew );
    tcdfDifference = log_tcdf_tauNew_nuNew - log_tcdf_tau_nu;

    log_y = - sum( log( abs( diag( cholOmegaCheck ) ) ) ) + MVTStudentTLogPDF_TIcholOmegaCheck_mInnovation_nu;

    FiniteCDFDifferenceSelect = isfinite( real( log_tcdf_tau_nu ) ) | isfinite( real( log_tcdf_tauNew_nuNew ) );
    
    log_y( FiniteCDFDifferenceSelect ) = log_y( FiniteCDFDifferenceSelect ) + tcdfDifference( FiniteCDFDifferenceSelect );
        
    assert( all( ~isnan( log_y(:) ) ), 'ESTNLSS:ESTLogPDF:NaNOutputLogY', 'ESTLogPDF returned a NaN output log_y.' );    
    
end
