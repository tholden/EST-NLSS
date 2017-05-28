function y = cbetainc( x, a, b )

    [ f, B ] = inbeta( x, a, b, 40 );
    y = f ./ B;

end
