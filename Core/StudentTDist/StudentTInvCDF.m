function x = StudentTInvCDF( y, nu )

    assert( numel( nu ) == 1, 'ESTNLSS:StudentTInvCDF:NuSize', 'StudentTInvCDF only supports univariate nu.' );
    assert( nu > 0, 'ESTNLSS:StudentTInvCDF:NuSign', 'StudentTInvCDF requires nu to be strictly positive.' );
    assert( all( ~isnan( y(:) ) ), 'ESTNLSS:StudentTInvCDF:NaNInputY', 'StudentTInvCDF was passed a NaN input y.' );
    
    x = tinv( y, nu );
    
    SelBad = ( x == -Inf ) & ( y > 0 );

    log_y = log( y );
    
    x( SelBad ) = -exp( 0.5 * log( nu ) - ( betaln( 0.5, 0.5 * nu ) + log( nu ) + log_y( SelBad ) ) / nu );
    
    % Use Newton's algorithm to polish
    Converged = ~isfinite( x );
    
    if ~all( Converged(:) )
    
        odx = Inf( size( x ) );
        for iter = 1 : 20
            
            [ ~, log_y_x ] = StudentTCDF( x, nu );
            [ ~, log_dy_x ] = StudentTPDF( x, nu );
            
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

    assert( all( ~isnan( x(:) ) ), 'ESTNLSS:StudentTInvCDF:NaNOutputX', 'StudentTInvCDF returned a NaN output x.' );
    
end
