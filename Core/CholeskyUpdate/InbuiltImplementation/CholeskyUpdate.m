function [ R, p ] = CholeskyUpdate( R, x, SignString )
    if nargin < 3
        SignString = '+';
    end
    if ( nargout < 2 ) || strcmp( SignString, '+' )
        p = 0;
        R = cholupdate( R, x, SignString );
    else
        [ R, p ] = cholupdate( R, x, SignString );
    end
end

% while true; R = randn( 10 ) + randn( 10 ) .* 1i; A = R' * R; R = chol( A ); R = R + diag( randn( 10, 1 ) .* 1i ); A = R' * R; x = randn( 10, 1 ) + randn( 10, 1 ) .* 1i; [ R1, p ] = CholeskyUpdate( R, x, '+' ); if p == 0; disp( max( max( abs( R1' * R1 - A - x * x' ) ) ) ); end; end;
% while true; R = randn( 10 ) + randn( 10 ) .* 1i; A = R' * R; R = chol( A ); R = R + diag( randn( 10, 1 ) .* 1i ); A = R' * R; x = randn( 10, 1 ) + randn( 10, 1 ) .* 1i; [ R1, p ] = CholeskyUpdate( R, x, '-' ); if p == 0; disp( max( max( abs( R1' * R1 - A + x * x' ) ) ) ); end; end;
