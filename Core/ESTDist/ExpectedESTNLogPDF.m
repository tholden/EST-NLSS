function [ log_y, Dlog_y ] = ExpectedESTNLogPDF( p, X, W, fInf, DynamicNu, SkewLikelihood, nu, mu, CholSigma )

    n = size( X, 1 );
    
    [ xi, CholOmega, delta, tau, nu ] = GetESTParametersFromVector( p, n, DynamicNu, SkewLikelihood, nu, mu, CholSigma );
    
    log_y_Points = -ESTLogPDF( X, xi, CholOmega, delta, tau, nu, true );
    log_y = sum( log_y_Points * W(:) );
    
    if ( ~isfinite( real( log_y ) ) ) || ( real( log_y ) > fInf )
        log_y = fInf;
    end
    
    if nargout > 1
        k = length( p );
        Dlog_y = zeros( 1, k );
        p = complex( p );
        se = sqrt( eps );
        sei = se * 1i;
        for i = 1 : k
            pTmp = p( i );
            p( i ) = pTmp + sei;
            Dlog_y( i ) = imag( ExpectedESTNLogPDF( p, X, W, fInf, DynamicNu, SkewLikelihood, nu, mu, CholSigma ) ) / se;
            p( i ) = pTmp;
        end
    end
    
end
