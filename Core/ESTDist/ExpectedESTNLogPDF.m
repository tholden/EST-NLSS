function [ log_y, Dlog_y ] = ExpectedESTNLogPDF( p, X, W )

    n = size( X, 1 );
    
    [ xi, CholOmega, delta, tau, nu ] = GetESTParametersFromVector( p, n );
    
    log_y_Points = -ESTLogPDF( X, xi, CholOmega, delta, tau, nu, true );
    
    log_y = log_y_Points * W(:);
    
    if nargout > 1
        k = length( p );
        Dlog_y = zeros( 1, k );
        p = complex( p );
        se = sqrt( eps );
        sei = se * 1i;
        for i = 1 : k
            pTmp = p( i );
            p( i ) = pTmp + sei;
            Dlog_y( i ) = imag( ExpectedESTNLogPDF( p, X, W ) ) / se;
            p( i ) = pTmp;
        end
    end
    
end
