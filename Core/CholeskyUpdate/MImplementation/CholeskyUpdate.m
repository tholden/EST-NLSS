function [ R, p ] = CholeskyUpdate( R, x, SignString )
    if nargin < 3 || strcmp( SignString, '+' )
        Sign = +1;
    else
        Sign = -1;
    end
    % derived from https://en.wikipedia.org/wiki/Cholesky_decomposition#Rank-one_update
    x = x(:);
    n = length( x );
    for p = 1 : n
        r2 = R( p, p ) .* R( p, p ) + Sign * x( p ) .* x( p );
        if r2 <= 0 || R( p, p ) == 0
            if nargout < 2
                error( 'ESTNLSS:NonPDFollowingCholeskyUpdate', 'CholeskyUpdate produced a non-positive-definite matrix.' );
            else
                return;
            end
        end
        r = sqrt( r2 );
        c = r / R( p, p );
        s = x( p ) / R( p, p );
        R( p, p ) = r;
        Indices = ( ( p + 1 ) : n )';
        R( p, Indices ) = ( R( p, Indices ) + Sign * s * x( Indices )' ) / c;
        x( Indices ) = c * x( Indices ) - s * R( p, Indices )';
    end
    p = 0;
end

% while true; R = randn( 10 ); A = R' * R; R = chol( A ); x = randn( 10, 1 ); [ R1, p ] = CholeskyUpdate( R, x, '+' ); if p == 0; disp( max( max( abs( R1' * R1 - A - x * x' ) ) ) ); end; end;
% while true; R = randn( 10 ); A = R' * R; R = chol( A ); x = randn( 10, 1 ); [ R1, p ] = CholeskyUpdate( R, x, '-' ); if p == 0; disp( max( max( abs( R1' * R1 - A + x * x' ) ) ) ); end; end;
