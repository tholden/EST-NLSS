function [ CubatureWeights, CubaturePoints, NCubaturePoints ] = GetGaussianCubaturePoints( IntDim, FilterCubatureDegree )
    if imag( FilterCubatureDegree ) ~= 0
        NCubaturePoints = 2 ^ imag( FilterCubatureDegree );
        CubaturePoints = [ zeros( IntDim, 1 ), randn( IntDim, NCubaturePoints - 1 ) ];
        CubatureWeights = ones( 1, NCubaturePoints );
    elseif FilterCubatureDegree < 0
        CubatureOrder = -FilterCubatureDegree;
        CubaturePoints = HigherOrderSobol( IntDim, CubatureOrder );
        NCubaturePoints = size( CubaturePoints, 2 );
        CubatureWeights = ones( 1, NCubaturePoints );
    elseif FilterCubatureDegree > 0
        CubatureOrder = ceil( 0.5 * ( FilterCubatureDegree - 1 ) );
        [ CubatureWeights, CubaturePoints, NCubaturePoints ] = fwtpts( IntDim, CubatureOrder );
    else
        NCubaturePoints = 2 * IntDim + 1;
        wTemp = 0.5 * realsqrt( 2 * NCubaturePoints );
        CubaturePoints = [ zeros( IntDim, 1 ), wTemp * eye( IntDim ), -wTemp * eye( IntDim ) ];
        CubatureWeights = ones( 1, NCubaturePoints );
    end
end
