addpath ../../Core/Utils

OrderMax = 20;
Dimension = 1;

CorrectMoments = [ sqrt( 2 / pi ), 1, 2 * sqrt( 2 / pi ), 3 ];

Moment = 4;

Correct = CorrectMoments( Moment );

LogErrors = NaN( OrderMax, 50 );

for Size = 1 : OrderMax
    for Smoothness = 1 : max( 1, min( 50, floor( 52 / Size ) ) )
        X = HigherOrderSobol( Dimension, Size, Smoothness );
        X = X( 1, : );
        LogErrors( Size, Smoothness ) = log( abs( mean( abs( X ) .^ Moment ) - Correct ) );
    end
end

figure( 1 ); plot( LogErrors, 'Marker', 'd' );
figure( 2 ); plot( LogErrors.', 'Marker', 'd' );
