function [ xi, InvCholOmega, delta, tau, nu, LogDetInvCholOmega ] = GetESTParametersFromVector( p, n )

    tau = p( 1 );
    nu = exp( p( 2 ) );
    xi = p( 3 : ( n + 2 ) );
    delta = p( ( n + 3 ) : ( 2 * n + 2 ) );
    LogDiagInvCholOmega = p( ( 2 * n + 3 ) : ( 3 * n + 2 ) );
    LTInvCholOmega = p( ( 3 * n + 3 ) : ( 3 * n + 2 + 0.5 * n * ( n - 1 ) ) );

    LogDetInvCholOmega = sum( LogDiagInvCholOmega );
    DiagInvCholOmega = exp( LogDiagInvCholOmega );
    
    InvCholOmega = tril( ones( n ), - 1 );
    InvCholOmega( ~~InvCholOmega ) = LTInvCholOmega;
    InvCholOmega = InvCholOmega + diag( DiagInvCholOmega );

end

