function x = SobolSetWrapper( Dim, n )
    ESTNLSSassert( n >= 1 + Dim, 'ESTNLSS:SobolSetWrapper:Dimensions', 'Input dimensions did not make sense.' );
    persistent p d
    if isempty( p ) || length( p ) < Dim || isempty( p{ Dim } )
        p{ Dim } = sobolset( Dim, 'Skip', 1 );
    end
    if isempty( d ) || length( d ) < Dim || isempty( d{ Dim } ) || size( d{ Dim }, 1 ) < n
        d{ Dim } = norminv( net( p{ Dim }, n ) )';
    end
    x = d{ Dim }( :, 1:n );
end
