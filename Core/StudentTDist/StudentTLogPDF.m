function log_y = StudentTLogPDF( x, nu )

    % To see why we do not use the MATLAB function, try plot( tpdf( 0, exp( 0:1:300 ) ) ) and plot( tpdf( 0, exp( 0:1:1000 ) ) )

    ESTNLSSassert( numel( nu ) == 1, 'ESTNLSS:StudentTLogPDF:NuSize', 'StudentTLogPDF only supports univariate nu.' );
    ESTNLSSassert( real( nu ) >= 0, 'ESTNLSS:StudentTLogPDF:NuSign', 'StudentTLogPDF requires nu to be weakly positive.' );
    ESTNLSSassert( all( ~isnan( x(:) ) ), 'ESTNLSS:StudentTLogPDF:NaNInputX', 'StudentTLogPDF was passed a NaN input x.' );

    if real( nu ) < Inf
        if real( nu ) > 0
            log_y = GetLogGammaRootNuRatio( nu ) - 0.6931471805599453094172 - 0.5 * ( nu + 1 ) .* log1p( ( x .* x ) ./ nu ); % -0.6931471805599453094172 = log(1/2)
        else
            log_y = -Inf;
        end
    else
        log_y = - 0.5 * ( x .* x ) - 0.918938533204672742;
    end
    
    ESTNLSSassert( all( ~isnan( log_y(:) ) ), 'ESTNLSS:StudentTLogPDF:NaNOutputLogY', 'StudentTLogPDF returned a NaN output log_y.' );    
    
end
