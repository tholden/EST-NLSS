function [ R, p ] = RealCholeskyUpdate( Rin, xin, SignString )
    if nargin < 3 || strcmp( SignString, '+' )
        Sign = +1;
    else
        Sign = -1;
    end
    % derived from https://en.wikipedia.org/wiki/Cholesky_decomposition#Rank-one_update
    xin = xin(:);
    n = length( xin );
    
    if isreal( xin ) || ( ~isreal( Rin ) )
        R = Rin;
    else
        R = complex( Rin );
    end
    
    if isreal( Rin ) || ( ~isreal( xin ) )
        x = xin;
    else
        x = complex( xin );
    end
    
    for k = 1 : n
        Rkk = R( k, k );
        xk = x( k );
        r2 = Rkk .* Rkk + Sign * xk .* xk;
        if real( r2 ) <= 0 || Rkk == 0
            if nargout < 2
                ESTNLSSerror( 'ESTNLSS:NonPDFollowingRealCholeskyUpdate', 'RealCholeskyUpdate produced a non-positive-definite matrix.' );
            else
                p = k;
                return
            end
        end
        r = sqrt( r2 );
        c = r / Rkk;
        s = xk / Rkk;
        R( k, k ) = r;
        Indices = ( ( k + 1 ) : n ).';
        R( k, Indices ) = ( R( k, Indices ) + Sign * s * x( Indices ).' ) / c;
        x( Indices ) = c * x( Indices ) - s * R( k, Indices ).';
    end
    p = 0;
end

% while true; R = randn( 10 ); A = R.' * R; R = chol( A ); x = randn( 10, 1 ); [ R1, p ] = RealCholeskyUpdate( R, x, '+' ); if p == 0; disp( max( max( abs( R1.' * R1 - A - x * x.' ) ) ) ); end; end;
% while true; R = randn( 10 ); A = R.' * R; R = chol( A ); x = randn( 10, 1 ); [ R1, p ] = RealCholeskyUpdate( R, x, '-' ); if p == 0; disp( max( max( abs( R1.' * R1 - A + x * x.' ) ) ) ); end; end;
