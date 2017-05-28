function Seq = Shanks( DSeq )

    Seq = cumsum( DSeq );
    
    if size( Seq, 1 ) >= 3
    
        DSeq = DSeq( 2:end, : );

        if mod( size( Seq, 1 ), 2 ) == 0
            Seq = Seq( 2:end, : );
            DSeq = DSeq( 2:end, : );
        end

        Seq = Seq( 3:end, : );
        LDSeq = DSeq( 1:(end-1), : );
        DSeq = DSeq( 2:end, : );

        Adj = DSeq .* DSeq ./ ( DSeq - LDSeq );
        FiniteAdj = isfinite( Adj );

        Seq( FiniteAdj ) = Seq( FiniteAdj ) - Adj( FiniteAdj );

        coder.unroll( );
        for Idx = 1 : ( 0.5 * ( size( Seq, 1 ) - 1 ) )

            DSeq = diff( Seq );

            Seq = Seq( 3:end, : );
            LDSeq = DSeq( 1:(end-1), : );
            DSeq = DSeq( 2:end, : );

            Adj = DSeq .* DSeq ./ ( DSeq - LDSeq );
            FiniteAdj = isfinite( Adj );

            Seq( FiniteAdj ) = Seq( FiniteAdj ) - Adj( FiniteAdj );

        end
    
    end
    
    Seq = Seq( end, : );
    
end
