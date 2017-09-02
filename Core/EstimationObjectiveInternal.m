function PersistentStateInternal0 = EstimationObjectiveInternal( StatDistPoints, PersistentStateInternal0, StateSteadyState, DynamicNu, SkewLikelihood, nuoo )
    if any( isnan( PersistentStateInternal0 ) )
    
        MeanStatDist = mean( StatDistPoints, 2 );
        DeMeanedStatDistPoints = bsxfun( @minus, StatDistPoints, MeanStatDist );
        [ ~, CholVarianceStatDist ] = NearestSPD( cov( StatDistPoints.' ) );

        MeanStatDistMMedianStatDist = MeanStatDist - StateSteadyState;
        CholVarianceStatDist_MeanStatDistMMedianStatDist = CholVarianceStatDist * MeanStatDistMMedianStatDist;
        CholVarianceStatDist_MeanStatDistMMedianStatDist2 = CholVarianceStatDist_MeanStatDistMMedianStatDist.' * CholVarianceStatDist_MeanStatDistMMedianStatDist;

        if CholVarianceStatDist_MeanStatDistMMedianStatDist2 > eps && SkewLikelihood
            ZcheckStatDist = ( MeanStatDistMMedianStatDist.' * DeMeanedStatDistPoints ) / realsqrt( CholVarianceStatDist_MeanStatDistMMedianStatDist2 );

            sZ3 = skewness( ZcheckStatDist, 0 );
            sZ4 = max( 3, kurtosis( ZcheckStatDist, 0 ) );

            if DynamicNu
                tauoo_nuoo = LMFnlsq2( @( in ) CalibrateMomentsEST( in( 1 ), 4 + eps( 4 ) + exp( in( 2 ) ), MeanStatDist, StateSteadyState, CholVarianceStatDist, sZ3, sZ4 ), [ 2; 2 ] );
                tauoo = tauoo_nuoo( 1 );
                nuoo = 4 + eps( 4 ) + exp( tauoo_nuoo( 2 ) );
            else
                tauoo = LMFnlsq2( @( in ) CalibrateMomentsEST( in( 1 ), nuoo, MeanStatDist, StateSteadyState, CholVarianceStatDist, sZ3, [] ), 2 );
            end
        else
            tauoo = Inf;

            if DynamicNu
                kurtDir = max( 0, kurtosis( DeMeanedStatDistPoints, 0, 2 ) - 3 );

                if kurtDir.' * kurtDir < eps
                    kurtDir = kurtosis( DeMeanedStatDistPoints, 0, 2 );
                end

                kurtDir = kurtDir / norm( kurtDir );

                ZcheckStatDist = kurtDir.' * DeMeanedStatDistPoints;

                sZ4 = max( 3, kurtosis( ZcheckStatDist, 0 ) );
                nuoo = 4 + 6 / ( sZ4 - 3 );
            end
        end

        [ ~, xoo, deltasoo, CholPsoo ] = CalibrateMomentsEST( tauoo, nuoo, MeanStatDist, StateSteadyState, CholVarianceStatDist, [], [] );
    
        PersistentStateInternal0 = InvGetESTParametersFromVector( xoo, CholPsoo, deltasoo, tauoo, nuoo, DynamicNu, SkewLikelihood );
        
    end
    
    W = 1 ./ size( StatDistPoints, 2 );
    
    p0 = PersistentStateInternal0;
    f0 = ExpectedESTNLogPDF( p0, StatDistPoints, W, Inf, DynamicNu, SkewLikelihood, nuoo );

    if coder.target( 'MATLAB' )
        fminlbfgsOptions = struct( 'Display', 'iter', 'GradObj', 'on', 'GradConstr', true, 'GoalsExactAchieve', false, 'TolX', 1e-12, 'TolFun', 1e-12, 'MaxIter', Inf, 'MaxFunEvals', Inf, ...
            'HessUpdate', 'bfgs ', 'DiffMaxChange', 1e-1, 'DiffMinChange', 1e-8, 'OutputFcn', [], 'rho', 0.0100, 'sigma', 0.900, 'tau1', 3, 'tau2', 0.1, 'tau3', 0.5, 'StoreN', 20 );
    else
        fminlbfgsOptions = struct( 'Display', 'off', 'GradObj', 'on', 'GradConstr', true, 'GoalsExactAchieve', false, 'TolX', 1e-12, 'TolFun', 1e-12, 'MaxIter', Inf, 'MaxFunEvals', Inf, ...
            'HessUpdate', 'bfgs ', 'DiffMaxChange', 1e-1, 'DiffMinChange', 1e-8, 'OutputFcn', [], 'rho', 0.0100, 'sigma', 0.900, 'tau1', 3, 'tau2', 0.1, 'tau3', 0.5, 'StoreN', 20 );
    end
    [ pOpt, fOpt ] = fminlbfgs( @( p ) ExpectedESTNLogPDF( p, StatDistPoints, W, f0 + 1, DynamicNu, SkewLikelihood, nuoo ), p0, fminlbfgsOptions );
    
    if fOpt < f0
        PersistentStateInternal0 = pOpt;
    else
        PersistentStateInternal0 = p0;
    end
end
