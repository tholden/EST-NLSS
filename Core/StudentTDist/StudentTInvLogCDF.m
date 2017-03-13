function x = StudentTInvLogCDF( log_y, nu )

    assert( numel( nu ) == 1, 'ESTNLSS:StudentTInvLogCDF:NuSize', 'StudentTInvLogCDF only supports univariate nu.' );
    assert( nu > 0, 'ESTNLSS:StudentTInvLogCDF:NuSign', 'StudentTInvLogCDF requires nu to be strictly positive.' );
    assert( all( ~isnan( log_y(:) ) ), 'ESTNLSS:StudentTInvLogCDF:NaNInputY', 'StudentTInvLogCDF was passed a NaN input log_y.' );
    
    x = tinv( exp( log_y ), nu );
    
    SelBad = ( x == -Inf ) & ( log_y > -Inf );

    x( SelBad ) = -exp( 0.5 * reallog( nu ) - ( betaln( 0.5, 0.5 * nu ) + reallog( nu ) + log_y( SelBad ) ) / nu );
    
    % Use Newton's algorithm to polish
    Converged = ~isfinite( x );
    
    if ~all( Converged(:) )
    
        odx = Inf( size( x ) );
        for iter = 1 : 20
            
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
    
    assert( all( ~isnan( x(:) ) ), 'ESTNLSS:StudentTInvLogCDF:NaNOutputX', 'StudentTInvLogCDF returned a NaN output x.' );
    
end
