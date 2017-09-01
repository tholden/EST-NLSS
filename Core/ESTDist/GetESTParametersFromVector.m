function [ xi, CholOmega, delta, tau, nu ] = GetESTParametersFromVector( p, n, DynamicNu, SkewLikelihood, nu )

    l = 2 * n + 0.5 * n * ( n - 1 );

    xi = p( 1 : n );
    LogDiagCholOmega = p( ( n + 1 ) : ( 2 * n ) );
    UpperCholOmega = p( ( 2 * n + 1 ) : l );

    if SkewLikelihood
        delta = p( ( l + 1 ) : ( l + n ) );
        tau = p( l + n + 1 );
    else
        delta = zeros( n, 1 );
        tau = Inf;
    end
    if DynamicNu
        nu = 2 + exp( p( end ) );
    end

    DiagCholOmega = exp( LogDiagCholOmega );
    
    if isreal( p )
        CholOmega = triu( ones( n ), 1 );
    else
        CholOmega = complex( triu( ones( n ), 1 ) );
    end
    CholOmega( CholOmega == 1 ) = UpperCholOmega;
    CholOmega = CholOmega + diag( DiagCholOmega );

end
