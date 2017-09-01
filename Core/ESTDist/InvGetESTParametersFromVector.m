function p = InvGetESTParametersFromVector( xi, CholOmega, delta, tau, nu, DynamicNu, SkewLikelihood )

    n = length( xi );
    l = 2 * n + 0.5 * n * ( n - 1 );
    p = zeros( l + ( n + 1 ) * SkewLikelihood + DynamicNu, 1 );
    
    LogDiagCholOmega = reallog( diag( CholOmega ) );
    UpperTriangle = triu( ones( n ), 1 );
    UpperCholOmega = CholOmega( UpperTriangle == 1 );
    
    p( 1 : n ) = xi;
    p( ( n + 1 ) : ( 2 * n ) ) = LogDiagCholOmega;
    p( ( 2 * n + 1 ) : l ) = UpperCholOmega;

    if SkewLikelihood
        p( ( l + 1 ) : ( l + n ) ) = delta;
        p( l + n + 1 ) = tau;
    end
    if DynamicNu
        p( end ) = log( nu - 2 );
    end

end
