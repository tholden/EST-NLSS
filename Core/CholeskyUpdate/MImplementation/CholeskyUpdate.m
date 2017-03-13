function [ R, p ] = CholeskyUpdate( R, x, SignString )
    if nargin < 3 || strcmp( SignString, '+' )
        Sign = +1;
    else
        Sign = -1;
    end
    % derived from https://en.wikipedia.org/wiki/Cholesky_decomposition#Rank-one_update
    x = x(:);
    n = length( x );
    for k = 1 : n
        r2 = R( k, k ) .* R( k, k ) + Sign * x( k ) .* x( k );
        if r2 <= 0 || R( k, k ) == 0
            if nargout < 2
                error( 'ESTNLSS:NonPDFollowingCholeskyUpdate', 'CholeskyUpdate produced a non-positive-definite matrix.' );
            else
                p = k;
                return;
            end
        end
        r = realsqrt( r2 );
        c = r / R( k, k );
        s = x( k ) / R( k, k );
        R( k, k ) = r;
        Indices = ( ( k + 1 ) : n )';
        R( k, Indices ) = ( R( k, Indices ) + Sign * s * x( Indices )' ) / c;
        x( Indices ) = c * x( Indices ) - s * R( k, Indices )';
    end
    p = 0;
end

% while true; R = randn( 10 ); A = R' * R; R = chol( A ); x = randn( 10, 1 ); [ R1, p ] = CholeskyUpdate( R, x, '+' ); if p == 0; disp( max( max( abs( R1' * R1 - A - x * x' ) ) ) ); end; end;
% while true; R = randn( 10 ); A = R' * R; R = chol( A ); x = randn( 10, 1 ); [ R1, p ] = CholeskyUpdate( R, x, '-' ); if p == 0; disp( max( max( abs( R1' * R1 - A + x * x' ) ) ) ); end; end;
