function x = StudentTInvCDF( y, nu )

    assert( numel( nu ) == 1 );
    assert( nu > 0 );
    
    x = tinv( y, nu );
    
    SelBad = ( x == -Inf ) & ( y > 0 );

    log_y = log( y );
    
    x( SelBad ) = -exp( 0.5 * log( nu ) - ( betaln( 0.5, 0.5 * nu ) + log( nu ) + log_y( SelBad ) ) / nu );
    
    % Use Newton's algorithm to polish
    Converged = ~isfinite( x );
    odx = Inf( size( x ) );
    for iter = 1 : 20
        [ ~, log_y_x ] = StudentTCDF( x, nu );
        [ ~, log_dy_x ] = StudentTPDF( x, nu );
        f_x = log_y_x - log_y;
        inv_df_x = exp( log_y_x - log_dy_x );
        dx = f_x .* inv_df_x;
        Converged = Converged | ( abs( dx ) <= max( eps, eps( x ) ) ) | ( dx > odx );
        odx = dx;
        x = x - dx;
        if all( Converged )
            break;
        end
    end
    
end
