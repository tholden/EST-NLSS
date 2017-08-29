function Seq = Shanks( DSeq )

    Seq = cumsum( DSeq );
    
    N = size( Seq, 1 );
    e = eps;
    se = sqrt( eps );
    
    if N >= 3
    
        DSeq = DSeq( 2:end, : );

        if mod( size( Seq, 1 ), 2 ) == 0
            Seq = Seq( 2:end, : );
            DSeq = DSeq( 2:end, : );
            N = N - 1;
        end

        Seq = Seq( 3:end, : );
        LDSeq = DSeq( 1:(end-1), : );
        DSeq = DSeq( 2:end, : );
        N = N - 2;

        Top = DSeq .* DSeq;
        Bottom = DSeq - LDSeq;
        
        aTop = abs( Top );
        aBottom = abs( Bottom );
        
        GoodAdj = ( aBottom > e ) & ( max( aTop, aBottom ) > se );

        Seq( GoodAdj ) = Seq( GoodAdj ) - Top( GoodAdj ) ./ Bottom( GoodAdj );

        coder.unroll( );
        for Idx = 1 : ( 0.5 * ( N - 1 ) )

            DSeq = diff( Seq );

            Seq = Seq( 3:end, : );
            LDSeq = DSeq( 1:(end-1), : );
            DSeq = DSeq( 2:end, : );

            Top = DSeq .* DSeq;
            Bottom = DSeq - LDSeq;

            aTop = abs( Top );
            aBottom = abs( Bottom );

            GoodAdj = ( aBottom > e ) & ( max( aTop, aBottom ) > se );

            Seq( GoodAdj ) = Seq( GoodAdj ) - Top( GoodAdj ) ./ Bottom( GoodAdj );

        end
    
    end
    
    Seq = Seq( end, : );
    
end
