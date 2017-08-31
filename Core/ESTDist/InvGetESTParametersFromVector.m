function p = InvGetESTParametersFromVector( xi, CholOmega, delta, tau, nu )

    n = length( xi );
    p = zeros( 3 * n + 2 + 0.5 * n * ( n - 1 ), 1 );
    
    LogDiagCholOmega = reallog( diag( CholOmega ) );
    UpperTriangle = triu( ones( n ), 1 );
    UpperCholOmega = CholOmega( UpperTriangle == 1 );
    
    p( 1 ) = tau;
    p( 2 ) = log( nu - 2 );
    p( 3 : ( n + 2 ) ) = xi;
    p( ( n + 3 ) : ( 2 * n + 2 ) ) = delta;
    p( ( 2 * n + 3 ) : ( 3 * n + 2 ) ) = LogDiagCholOmega;
    p( ( 3 * n + 3 ) : ( 3 * n + 2 + 0.5 * n * ( n - 1 ) ) ) = UpperCholOmega;

end
