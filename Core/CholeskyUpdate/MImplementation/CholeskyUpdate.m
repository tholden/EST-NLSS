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
        Rkk = R( k, k );
        conjRkk = conj( Rkk );
        xk = x( k );
        r2 = conjRkk .* Rkk + Sign * conj( xk ) .* xk;
        if real( r2 ) <= 0 || Rkk == 0
            if nargout < 2
                ESTNLSSerror( 'ESTNLSS:NonPDFollowingCholeskyUpdate', 'CholeskyUpdate produced a non-positive-definite matrix.' );
            else
                p = k;
                return
            end
        end
        r = sqrt( r2 );
        c = r / conjRkk;
        s = xk / conjRkk;
        R( k, k ) = r;
        Indices = ( ( k + 1 ) : n )';
        R( k, Indices ) = ( R( k, Indices ) + Sign * s * x( Indices )' ) / c;
        x( Indices ) = c * x( Indices ) - s * R( k, Indices )';
    end
    p = 0;
end

% while true; R = randn( 10 ) + randn( 10 ) .* 1i; A = R' * R; R = chol( A ); R = R + diag( randn( 10, 1 ) .* 1i ); A = R' * R; x = randn( 10, 1 ) + randn( 10, 1 ) .* 1i; [ R1, p ] = CholeskyUpdate( R, x, '+' ); if p == 0; disp( max( max( abs( R1' * R1 - A - x * x' ) ) ) ); end; end;
% while true; R = randn( 10 ) + randn( 10 ) .* 1i; A = R' * R; R = chol( A ); R = R + diag( randn( 10, 1 ) .* 1i ); A = R' * R; x = randn( 10, 1 ) + randn( 10, 1 ) .* 1i; [ R1, p ] = CholeskyUpdate( R, x, '-' ); if p == 0; disp( max( max( abs( R1' * R1 - A + x * x' ) ) ) ); end; end;
