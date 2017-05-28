function y = cbetaln( a, b )

    y = cgammaln( a ) + cgammaln( b ) - cgammaln( a + b );

end
