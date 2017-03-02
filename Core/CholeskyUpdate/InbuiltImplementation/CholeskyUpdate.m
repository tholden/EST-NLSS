function [ R, p ] = CholeskyUpdate( R, x, SignString )
    if nargin < 3
        SignString = '+';
    end
    if nargout < 2
        R = cholupdate( R, x, SignString );
    else
        [ R, p ] = cholupdate( R, x, SignString );
    end
end

% while true; R = randn( 10 ); A = R' * R; R = chol( A ); x = randn( 10, 1 ); [ R1, p ] = CholeskyUpdate( R, x, '+' ); if p == 0; disp( max( max( abs( R1' * R1 - A - x * x' ) ) ) ); end; end;
% while true; R = randn( 10 ); A = R' * R; R = chol( A ); x = randn( 10, 1 ); [ R1, p ] = CholeskyUpdate( R, x, '-' ); if p == 0; disp( max( max( abs( R1' * R1 - A + x * x' ) ) ) ); end; end;
