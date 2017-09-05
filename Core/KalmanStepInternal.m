function PersistentStateInternal_t = KalmanStepInternal( wmRed, CubatureWeights, PersistentStateInternal_t, DynamicNu, SkewLikelihood, nuno, mu_wmRed, CholSigma_wmRed )

    if any( isnan( PersistentStateInternal_t ) )
        % Then initialize by calibration
        
        Median_wmRed = wmRed( :, 1 );

        Demeaned_wmRed = bsxfun( @minus, wmRed, mu_wmRed );
        
        Mean_wmRedMMedian_wmRed = mu_wmRed - Median_wmRed;
        CholVariance_wmRed_Mean_wmRedMMedian_wmRed = CholSigma_wmRed * Mean_wmRedMMedian_wmRed;
        CholVariance_wmRed_Mean_wmRedMMedian_wmRed2 = CholVariance_wmRed_Mean_wmRedMMedian_wmRed.' * CholVariance_wmRed_Mean_wmRedMMedian_wmRed;

        if CholVariance_wmRed_Mean_wmRedMMedian_wmRed2 > eps && SkewLikelihood
            Zcheck_wmRed = ( Mean_wmRedMMedian_wmRed' * Demeaned_wmRed ) / realsqrt( CholVariance_wmRed_Mean_wmRedMMedian_wmRed2 );

            meanZcheck_wmRed = Zcheck_wmRed * CubatureWeights.';
            Zcheck_wmRed = Zcheck_wmRed - meanZcheck_wmRed;
            meanZcheck_wmRed2 = ( Zcheck_wmRed .* Zcheck_wmRed ) * CubatureWeights.';
            Zcheck_wmRed = Zcheck_wmRed / realsqrt( meanZcheck_wmRed2 );

            sZ3 = realpow( Zcheck_wmRed, 3 ) * CubatureWeights.';
            sZ4 = max( 3, realpow( Zcheck_wmRed, 4 ) * CubatureWeights.' );

            if DynamicNu
                tauno_nuno = LMFnlsq2( @( in ) CalibrateMomentsEST( in( 1 ), 5 + exp( in( 2 ) ), mu_wmRed, Median_wmRed, CholSigma_wmRed, sZ3, sZ4 ), [ 2; 0 ] );
                tauno = tauno_nuno( 1 );
                nuno = 5 + exp( tauno_nuno( 2 ) );
            else
                tauno = LMFnlsq2( @( in ) CalibrateMomentsEST( in( 1 ), nuno, mu_wmRed, Median_wmRed, CholSigma_wmRed, sZ3, [] ), 2 );
            end
        else
            tauno = Inf;

            if DynamicNu
                Zcheck_wmRed = Demeaned_wmRed;

                meanZcheck_wmRed2 = ( Zcheck_wmRed .* Zcheck_wmRed ) * CubatureWeights.';
                Zcheck_wmRed = bsxfun( @times, Zcheck_wmRed, 1 ./ realsqrt( meanZcheck_wmRed2 ) );

                kurtDir = max( 0, realpow( Zcheck_wmRed, 4 ) * CubatureWeights.' - 3 );

                if kurtDir.' * kurtDir < eps
                    kurtDir = realpow( Zcheck_wmRed, 4 ) * CubatureWeights.';
                end

                kurtDir = kurtDir / norm( kurtDir );

                Zcheck_wmRed = kurtDir.' * Zcheck_wmRed;

                meanZcheck_wmRed = Zcheck_wmRed * CubatureWeights.';
                Zcheck_wmRed = Zcheck_wmRed - meanZcheck_wmRed;
                meanZcheck_wmRed2 = ( Zcheck_wmRed .* Zcheck_wmRed ) * CubatureWeights.';
                Zcheck_wmRed = Zcheck_wmRed / realsqrt( meanZcheck_wmRed2 );

                sZ4 = max( 3, realpow( Zcheck_wmRed, 4 ) * CubatureWeights.' );
                nuno = 4 + 6 / ( sZ4( 1 ) - 3 );
            end
        end

        if SkewLikelihood
            tauno = max( -10, min( 10, tauno ) );
        end
        if DynamicNu
            nuno = max( 5 + sqrt( eps( 5 ) ), min( 1000, nuno ) );
        end
        
        [ ~, ~, deltaetano ] = CalibrateMomentsEST( tauno, nuno, mu_wmRed, Median_wmRed, CholSigma_wmRed, [], [] );
        
        PersistentStateInternal_t = InvGetESTParametersFromVector( deltaetano, tauno, nuno, DynamicNu, SkewLikelihood );
    end
    
    p0 = PersistentStateInternal_t;
    f0 = ExpectedESTNLogPDF( p0, wmRed, CubatureWeights, Inf, DynamicNu, SkewLikelihood, nuno, mu_wmRed, CholSigma_wmRed );
    
    if coder.target( 'MATLAB' )
        fminlbfgsOptions = struct( 'Display', 'final', 'GradObj', 'on', 'GradConstr', true, 'GoalsExactAchieve', false, 'TolX', 1e-12, 'TolFun', 1e-12, 'MaxIter', 200 + 10 * ( length( p0 ) .^ 2 ), 'MaxFunEvals', Inf, ...
            'HessUpdate', 'bfgs ', 'DiffMaxChange', 1e-1, 'DiffMinChange', 1e-8, 'OutputFcn', [], 'rho', 0.0100, 'sigma', 0.900, 'tau1', 3, 'tau2', 0.1, 'tau3', 0.5, 'StoreN', 20 );
    else
        fminlbfgsOptions = struct( 'Display', 'off', 'GradObj', 'on', 'GradConstr', true, 'GoalsExactAchieve', false, 'TolX', 1e-12, 'TolFun', 1e-12, 'MaxIter', 200 + 10 * ( length( p0 ) .^ 2 ), 'MaxFunEvals', Inf, ...
            'HessUpdate', 'bfgs ', 'DiffMaxChange', 1e-1, 'DiffMinChange', 1e-8, 'OutputFcn', [], 'rho', 0.0100, 'sigma', 0.900, 'tau1', 3, 'tau2', 0.1, 'tau3', 0.5, 'StoreN', 20 );
    end
    [ pOpt, fOpt ] = fminlbfgs( @( p ) ExpectedESTNLogPDF( p, wmRed, CubatureWeights, f0 + 1, DynamicNu, SkewLikelihood, nuno, mu_wmRed, CholSigma_wmRed ), p0, fminlbfgsOptions );
    if fOpt < f0
        PersistentStateInternal_t = pOpt;
    else
        PersistentStateInternal_t = p0;
    end
end
