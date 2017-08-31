function [ xi, CholOmega, delta, tau, nu ] = GetESTParametersFromVector( p, n )

    tau = p( 1 );
    nu = 2 + exp( p( 2 ) );
    xi = p( 3 : ( n + 2 ) );
    delta = p( ( n + 3 ) : ( 2 * n + 2 ) );
    LogDiagCholOmega = p( ( 2 * n + 3 ) : ( 3 * n + 2 ) );
    UpperCholOmega = p( ( 3 * n + 3 ) : ( 3 * n + 2 + 0.5 * n * ( n - 1 ) ) );

    DiagCholOmega = exp( LogDiagCholOmega );
    
    CholOmega = triu( ones( n ), 1 );
    CholOmega( CholOmega == 1 ) = UpperCholOmega;
    CholOmega = CholOmega + diag( DiagCholOmega );

end
