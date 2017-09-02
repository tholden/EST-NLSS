function log_y = MVTStudentTLogPDF( x, nu )

    ESTNLSSassert( numel( nu ) == 1, 'ESTNLSS:MVTStudentTLogPDF:NuSize', 'MVTStudentTLogPDF only supports univariate nu.' );
    ESTNLSSassert( real( nu ) > 0, 'ESTNLSS:MVTStudentTLogPDF:NuSign', 'MVTStudentTLogPDF requires nu to be strictly positive.' );
    ESTNLSSassert( all( ~isnan( x(:) ) ), 'ESTNLSS:MVTStudentTLogPDF:NaNInputX', 'MVTStudentTLogPDF was passed a NaN input x.' );
    
    D = size( x, 1 );
    HD = 0.5 * D;
    fHD = floor( HD );
    
    if isfinite( real( nu ) )
        nuPDO2 = 0.5 * ( nu + D );
        logGammaRatio = sum( log( 0.5 * nu + ( 0 : ( fHD - 1 ) ) ) );
        if HD ~= fHD
            fnu = nu + D - 1;
            if real( fnu ) > 20
                ifnu = 1 ./ fnu;
                ifnu2 = ifnu .* ifnu;
                ifnu4 = ifnu2 .* ifnu2;
                ifnu6 = ifnu4 .* ifnu2;
                ifnu8 = ifnu4 .* ifnu4;
                ifnu10 = ifnu6 .* ifnu4;
                ifnu12 = ifnu6 .* ifnu6;
                ifnu14 = ifnu8 .* ifnu6;
                ifnu16 = ifnu8 .* ifnu8;
                logGammaRatio = logGammaRatio + 0.5 * log( fnu ) -.346573590279972655 - ifnu .* ( ...
                    .250000000000000000 + 0.416666666666666667e-1.*ifnu2 - 0.500000000000000000e-1.*ifnu4 + .151785714285714286.*ifnu6 - .861111111111111111.*ifnu8 + 7.85227272727272727.*ifnu10 - 105.019230769230769.*ifnu12 + 1936.60208333333333.*ifnu14 - 47092.5147058823529.*ifnu16 ...
                );
            else
                if isreal( nu )
                    logGammaRatio = logGammaRatio + gammaln( nuPDO2 ) - gammaln( 0.5 * fnu );
                else
                    logGammaRatio = logGammaRatio + cgammaln( nuPDO2 ) - cgammaln( 0.5 * fnu );
                end
            end
        end
        if isreal( nu )
            log_y = logGammaRatio - HD * reallog( nu ) - HD * 1.14472988584940017 - nuPDO2 .* log1p( sum( x .* x, 1 ) ./ nu ); % log( pi ) = 1.14472988584940017
        else
            log_y = logGammaRatio - HD * log( nu ) - HD * 1.14472988584940017 - nuPDO2 .* log1p( sum( x .* x, 1 ) ./ nu ); % log( pi ) = 1.14472988584940017
        end
    else
        log_y = - 0.5 * sum( x .* x, 1 ) - D * 0.91893853320467274; % 0.5 * log( 2 * pi ) = 0.91893853320467274
    end
    
    ESTNLSSassert( all( ~isnan( log_y(:) ) ), 'ESTNLSS:MVTStudentTLogPDF:NaNOutputLogY', 'MVTStudentTLogPDF returned a NaN output log_y.' );    
    
end
