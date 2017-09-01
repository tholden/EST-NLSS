function Seq = Shanks( DSeq )

    Seq = cumsum( DSeq );
    
    [ N, D ] = size( Seq );
    e = eps;
    se = sqrt( eps );
    
    if N >= 3
    
        if mod( size( Seq, 1 ), 2 ) == 0
            SeqIdx = 2;
            DSeqIdx = 3;
            N = N - 1;
        else
            SeqIdx = 1;
            DSeqIdx = 2;
        end

        SeqIdx = SeqIdx + 2;
        DSeqIdx = DSeqIdx + 1;
        N = N - 2;

        if isreal( DSeq )
            Top = coder.nullcopy( zeros( size( Seq ) ) );
            Bottom = coder.nullcopy( zeros( size( Seq ) ) );
        else
            Top = coder.nullcopy( complex( zeros( size( Seq ) ) ) );
            Bottom = coder.nullcopy( complex( zeros( size( Seq ) ) ) );
        end
        
        Top( SeqIdx : end, : ) = DSeq( DSeqIdx : end, : ) .* DSeq( DSeqIdx : end, : );
        Bottom( SeqIdx : end, : ) = DSeq( DSeqIdx : end, : ) - DSeq( ( DSeqIdx - 1 ) : ( end - 1 ), : );
        
        GoodAdj = [ false( SeqIdx - 1, D ); ( abs( Bottom( SeqIdx : end, : ) ) > e ) & ( max( abs( Top( SeqIdx : end, : ) ), abs( Bottom( SeqIdx : end, : ) ) ) > se ) ];

        Seq( GoodAdj ) = Seq( GoodAdj ) - Top( GoodAdj ) ./ Bottom( GoodAdj );

        for Idx = coder.unroll( 1 : ( 0.5 * ( N - 1 ) ) )

            DSeqIdx = DSeqIdx + 1;
            DSeq( DSeqIdx : end, : ) = Seq( ( SeqIdx + 1 ) : end, : ) - Seq( SeqIdx : ( end - 1 ), : );

            SeqIdx = SeqIdx + 2;
            DSeqIdx = DSeqIdx + 1;

            Top( SeqIdx : end, : ) = DSeq( DSeqIdx : end, : ) .* DSeq( DSeqIdx : end, : );
            Bottom( SeqIdx : end, : ) = DSeq( DSeqIdx : end, : ) - DSeq( ( DSeqIdx - 1 ) : ( end - 1 ), : );

            GoodAdj = [ false( SeqIdx - 1, D ); ( abs( Bottom( SeqIdx : end, : ) ) > e ) & ( max( abs( Top( SeqIdx : end, : ) ), abs( Bottom( SeqIdx : end, : ) ) ) > se ) ];

            Seq( GoodAdj ) = Seq( GoodAdj ) - Top( GoodAdj ) ./ Bottom( GoodAdj );

        end
    
    end
    
    Seq = Seq( end, : );
    
end
