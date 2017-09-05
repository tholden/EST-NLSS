function p = InvGetESTParametersFromVector( delta, tau, nu, DynamicNu, SkewLikelihood )

    n = length( delta );
    p = zeros( ( n + 1 ) * SkewLikelihood + DynamicNu, 1 );
    
    if SkewLikelihood
        p( 1 : n ) = delta;
        p( n + 1 ) = tau;
    end
    if DynamicNu
        p( end ) = log( nu - 5 );
    end

end
