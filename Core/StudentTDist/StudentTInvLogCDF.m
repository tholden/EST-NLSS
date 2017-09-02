function x = StudentTInvLogCDF( log_y, nu )

    ESTNLSSassert( numel( nu ) == 1, 'ESTNLSS:StudentTInvLogCDF:NuSize', 'StudentTInvLogCDF only supports univariate nu.' );
    ESTNLSSassert( nu > 0, 'ESTNLSS:StudentTInvLogCDF:NuSign', 'StudentTInvLogCDF requires nu to be strictly positive.' );
    ESTNLSSassert( all( ~isnan( log_y(:) ) ), 'ESTNLSS:StudentTInvLogCDF:NaNInputY', 'StudentTInvLogCDF was passed a NaN input log_y.' );
    
    x = tinv( exp( log_y ), nu );
    
    SelBad = ( x == -Inf ) & ( log_y > -Inf );

    betalnTmp = NaN;
    if isfinite( nu )
        betalnTmp = betaln( 0.5, 0.5 * nu );
        if isfinite( betalnTmp )
            x( SelBad ) = -exp( 0.5 * reallog( nu ) - ( betalnTmp + reallog( nu ) + log_y( SelBad ) ) / nu );
        end
    end
    if ~isfinite( betalnTmp )
        log_y_SelBad = log_y( SelBad );
        log_Mx_SelBad = zeros( size( log_y_SelBad ) );
        NormConst = 0.918938533204672741780329736407;
        for norm_iter = coder.unroll( 1 : 5 )
            Mx_SelBad = realsqrt( -2 * min( 0, log_y_SelBad + NormConst + log_Mx_SelBad ) );
            log_Mx_SelBad = log( Mx_SelBad );
        end
        x( SelBad ) = -Mx_SelBad;
    end
    
    % Use Newton's algorithm to polish
    Converged = ~isfinite( x );
    
    if ~all( Converged(:) )
    
        odx = Inf( size( x ) );
        for iter = coder.unroll( 1 : 20 )
            
            log_y_x = StudentTLogCDF( x, nu );
            log_dy_x = StudentTLogPDF( x, nu );
            
            f_x = log_y_x - log_y;
            inv_df_x = exp( log_y_x - log_dy_x );
            dx = f_x .* inv_df_x;
            
            Converged = Converged | ( abs( dx ) <= max( eps, eps( x ) ) ) | ( dx > odx ) | ~isfinite( dx );
            
            dx( ~isfinite( dx ) ) = 0;
            odx = dx;
            x = x - dx;
            
            if all( Converged(:) )
                break;
            end
            
        end
    
    end
    
    ESTNLSSassert( all( ~isnan( x(:) ) ), 'ESTNLSS:StudentTInvLogCDF:NaNOutputX', 'StudentTInvLogCDF returned a NaN output x.' );
    
end
