function X = HigherOrderSobol( Size, Dimension, Smoothness, AddZeroPoint )
    % Higher order Sobol sequence
    % Derived from https://quasirandomideas.wordpress.com/2010/06/17/how-to-generate-higher-order-sobol-points-in-matlab-and-some-numerical-examples/
    % Create a higher order Sobol sequence.
    % 2^Size:       number of points
    % Dimension:    dimension of final point set
    % Smoothness:   interlacing factor
    % AddZeroPoint: whether to add a point at 0
    % X Output Sobol sequence
    
    SmoothnessMax = min( 50, floor( 52 / Size ) );
    
    if nargin < 3
        Smoothness = SmoothnessMax;
    else
        Smoothness = min( Smoothness, SmoothnessMax );
    end
    
    if nargin < 4
        AddZeroPoint = false;
    end
    
    Smoothness = max( 1, Smoothness );

    persistent HigherOrderSobolCache

    % coder.varsize( 'X', [], [ true, true ] );
    % X = zeros( Dimension, 0 ); %#ok<PREALL>

    FileName = 'HigherOrderSobolCache.mat';
    
    if isempty( HigherOrderSobolCache )
        FoundFile = false;

        if exist( FileName, 'file' ) == 2
            LoadedFile = load( FileName );
            if isstruct( LoadedFile ) && isfield( LoadedFile, 'HigherOrderSobolCache' ) && ~isempty( LoadedFile.HigherOrderSobolCache )
                HigherOrderSobolCache = LoadedFile.HigherOrderSobolCache;
                FoundFile = true;
            end
        end
        
        if ~FoundFile
            % coder.varsize( 'HigherOrderSobolCache', [], [ true, false ] );
            HigherOrderSobolCache = struct( 'InputParameters', [ 0, 0, 0, 0 ], 'X', zeros( 0, 0 ) );
            % coder.varsize( 'HigherOrderSobolCache(:).X', [], [ true, true ] );
        end
    end

    InputParameters = [ Size, Dimension, Smoothness, AddZeroPoint ];
    for i = 1 : numel( HigherOrderSobolCache )
        if all( HigherOrderSobolCache( i ).InputParameters == InputParameters )
            X = HigherOrderSobolCache( i ).X;
            return
        end
    end
        
    N = pow2( Size ); % Number of points;
    P = sobolset( Smoothness * Dimension ); % Get Sobol sequence;
    SobolPoints = net( P, N ); % Get net from Sobol sequence with N points;

    % Create binary representation of digits;

    if Smoothness == 1
        X = SobolPoints.';
    else
        Z = SobolPoints.' * N;
        X = zeros( Dimension, N );
        for j = 1 : Dimension
            for i = 1 : Size
                for k = 1 : Smoothness
                    X( j, : ) = bitset( X(j,:), (Size*Smoothness+1) - k - (i-1)*Smoothness, bitget( Z((j-1)*Smoothness+k,:), (Size+1) - i ) );
                end
            end
        end
        X = X * pow2( -Size * Smoothness );
    end

    for i = 1 : Dimension

        Xi = X( i, : );

        SortedXi = sort( Xi );

        if Smoothness == 1
            NumPossible = 1;
        else
            NumPossible = pow2( floor( 0.5 * Size ) );
        end
        Possible = 1 : NumPossible;

        fprintf( '\n%d %d %d 1\n', int32( Size ), int32( Smoothness ), i );

        [ BestXCandidate, BestMinPenalty, BestLocation ] = PerformSearch( Xi, SortedXi, Possible, AddZeroPoint );

        if Smoothness > 1

            Possible = ( BestLocation + NumPossible ) : NumPossible : N;

            fprintf( '\n%d %d %d 2\n', int32( Size ), int32( Smoothness ), i );

            [ NewBestXCandidate, NewBestMinPenalty ] = PerformSearch( Xi, SortedXi, Possible, AddZeroPoint );

            if NewBestMinPenalty < BestMinPenalty
                BestXCandidate = NewBestXCandidate;
            end

        end

        X( i, : ) = BestXCandidate;

    end

    if AddZeroPoint
        X = [ zeros( Dimension, 1 ), X ];
    end

    assert( all( abs( mean( X, 2 ) ) < sqrt( eps ) ), 'ESTNLSS:HigherOrderSobol:Uncentered', 'Result was not centred.' );
    
    ToCache = struct( 'InputParameters', InputParameters, 'X', X );

    HigherOrderSobolCache = [ HigherOrderSobolCache; ToCache ];
    
    save( FileName, 'HigherOrderSobolCache' );

end

function [ BestXCandidate, BestMinPenalty, BestLocation ] = PerformSearch( Xi, SortedXi, Possible, AddZeroPoint )

    BatchSize = 64;
    
    NumPossible = numel( Possible );

    BestMinPenalty = Inf;
    BestXCandidate = [];
    BestLocation = 0;
    
    jMax = floor( ( NumPossible - 1 ) / BatchSize );

    for j = 0 : jMax

        fprintf( '%d / %d\n', int32( j + 1 ), int32( jMax + 1 ) );

        Current = j * BatchSize + ( 1 : min( NumPossible - j * BatchSize, BatchSize ) );
        NumCurrent = numel( Current );

        XiMat = bsxfun( @minus, Xi, SortedXi( Possible( Current ) )' );
        XiMat = XiMat - floor( XiMat );

        uiMax = 1 - max( XiMat, [], 2 );
        f = Inf( NumCurrent, 1 );
        vi = 0.5 * ones( NumCurrent, 1 );
        idxStillGoing = 1 : NumCurrent;

        k = 0;
        while ~isempty( idxStillGoing )
            of = f;
            vic = vi( idxStillGoing );
            [ f, dfdv, df2dv2 ] = GetResid( vic, XiMat( idxStillGoing, : ), uiMax( idxStillGoing ) );
            Step = max( -0.5 * vic, min( 0.5 * ( 1 - vic ), ( ( -2 ) .* ( f .* dfdv ) ) ./ ( 2 .* ( dfdv .* dfdv ) - f .* df2dv2 ) ) );
            vi( idxStillGoing ) = vic + Step;
            RelStillGoing = ( abs( Step ) > max( 1, abs( vi( idxStillGoing ) ) ) * eps ) & ( abs( f ) < abs( of ) );
            f = f( RelStillGoing );
            idxStillGoing = idxStillGoing( RelStillGoing );
            k = k + 1;
        end

        ui = uiMax .* vi;

        assert( all( isfinite( ui ) ), 'ESTNLSS:HigherOrderSobol:SolutionFailure', 'Failed to solve for a rotation.' );

        NewXCandidates = norminv( bsxfun( @plus, XiMat, ui ) );

        if AddZeroPoint
            TmpCandidates = [ zeros( NumCurrent, 1 ), NewXCandidates ];
        else
            TmpCandidates = NewXCandidates;
        end

        AbsSkewness = abs( mean( TmpCandidates .^ 3, 2 ) );
        AbsStdDevErr = abs( mean( TmpCandidates .^ 2, 2 ) - 1 );

        Penalty = 10000 * AbsSkewness + AbsStdDevErr;

        Penalty( ~isfinite( Penalty ) ) = max( Penalty( isfinite( Penalty ) ) );

        [ MinPenalty, IdxMinPenalty ] = min( Penalty );

        if MinPenalty < BestMinPenalty
            BestXCandidate = NewXCandidates( IdxMinPenalty, : );
            BestMinPenalty = MinPenalty;
            BestLocation = Possible( Current( IdxMinPenalty ) );
        end

    end

end

function [ f, dfdv, df2dv2 ] = GetResid( v, X, uMax )

    % N = size( X, 2 );

    u = uMax .* v;
    
    X = max( eps, min( 1 - eps, bsxfun( @plus, X, u ) ) );
    % dXdv = repmat( uMax, 1, N );
    % dX2dv2 = repmat( uMax, 1, N );
    
    X = -1.4142135623730950488 .* erfcinv( 2 * X ); % = norminv( X );
    tmp = exp( 0.5 .* ( X .* X ) );
    dXdv = 2.5066282746310005024 .* bsxfun( @times, uMax, tmp ); % = dxdv ./ normpdf( Y );
    dX2dv2 = dXdv .* ( 1 + X .* dXdv );
    
    f = sum( X, 2 );
    dfdv = sum( dXdv, 2 );
    df2dv2 = sum( dX2dv2, 2 );
    
end

%         ToTest = log2( NumPossible ) - 1;
%         AICcVec = zeros( ToTest, 1 );
%         BICcVec = zeros( ToTest, 1 );
%         MSEVec = zeros( ToTest, 1 );
%         assert( all( isfinite( LogPenaltyAll ) ) );
%         for j = 0 : ToTest
%             k = pow2( j );
%             TmpLogPenaltyAll = reshape( LogPenaltyAll(:), k, NumPossible / k );
%             TmpLogPenaltyAll1 = bsxfun( @minus, TmpLogPenaltyAll, mean( TmpLogPenaltyAll, 2 ) );
%             TmpLogPenaltyAll1 = TmpLogPenaltyAll1(:);
%             TmpLogPenaltyAll2 = bsxfun( @minus, TmpLogPenaltyAll( :, 2:end ), TmpLogPenaltyAll( :, 1 ) );
%             TmpLogPenaltyAll2 = TmpLogPenaltyAll2(:);
%             TmpRSS = sum( TmpLogPenaltyAll1 .* TmpLogPenaltyAll1 );
%             k = k + 1;
%             AICcVec( j + 1 ) = 2 * k + NumPossible * log( TmpRSS ) + 2 * k * ( k + 1 ) / ( NumPossible - k - 1 );
%             BICcVec( j + 1 ) = k * ( log( NumPossible ) - log( 2 * pi ) ) + NumPossible * log( TmpRSS );
%             MSEVec( j + 1 ) = mean( TmpLogPenaltyAll2 .* TmpLogPenaltyAll2 );
%         end
%         [ ~, AICcMinLoc ] = min( AICcVec );
%         [ ~, BICcMinLoc ] = min( BICcVec );
%         [ ~, MSEMinLoc ] = min( MSEVec );
%         fprintf( '%d %d %d %d %d %d\n', int32( Size ), int32( Smoothness ), int32( i ), int32( AICcMinLoc - 1 ), int32( BICcMinLoc - 1 ), int32( MSEMinLoc - 1 ) );
%         
%         figure( 3 ); plot( reshape( LogPenaltyAll(:), pow2( floor( 0.5 * log2( NumPossible ) ) ), NumPossible / pow2( floor( 0.5 * log2( NumPossible ) ) ) ) ); drawnow;
