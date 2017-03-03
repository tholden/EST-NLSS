function [ x, f, PersistentState ] = ACDWrapper( OptiFunction, x, lb, ub, PersistentState, varargin )

    try
        pool = gcp;
        nw = pool.NumWorkers;
    catch
        try
            nw = matlabpool( 'size' ); %#ok<DPOOL>
        catch
            nw = 1;
        end
    end

    Radius = 0.3;
    sigma = ( ub - lb ) * Radius;
    sigma( ~isfinite( sigma ) ) = max( Radius, Radius * abs( x( ~isfinite( sigma ) ) ) );
    
    InitialTimeOutLikelihoodEvaluation = Inf;
    
    ProductSearchDimension = max( 1, floor( sqrt( nw ) ) );
    Order = max( 1, ceil( ( -log( 2 ) + log( ( 20 * nw + 1 ) ^ ( 1 / ProductSearchDimension ) + 1 ) ) / log( 2 ) ) );
    
    [ x, f, PersistentState ] = ACDMinimisation( ...
        @( XV, PS, DesiredNumberOfNonTimeouts ) ParallelWrapper( @( X ) OptiFunction( X, PS ), XV, DesiredNumberOfNonTimeouts, InitialTimeOutLikelihoodEvaluation ),...
        x, sigma, eps, lb, ub, [], [], Inf, -Inf, 1, Order, 1, ProductSearchDimension, PersistentState, false );
    
    x = max( lb, min( ub, x ) );
    
    f = -f;
    
end

