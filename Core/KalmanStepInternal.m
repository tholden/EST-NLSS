function PersistentStateInternal_t = KalmanStepInternal( wmRed, CubatureWeights, PersistentStateInternal_t, DynamicNu, SkewLikelihood, nuno )

    if any( isnan( PersistentStateInternal_t ) )
        % Then initialize by calibration
        
        Median_wmRed = wmRed( :, 1 );

        Mean_wmRed = sum( bsxfun( @times, wmRed, CubatureWeights ), 2 );
        ano = bsxfun( @minus, wmRed, Mean_wmRed );
        Weighted_ano = bsxfun( @times, ano, CubatureWeights );

        Variance_wmRed = NearestSPD( ano * Weighted_ano.' );
        Variance_wmRed = 0.5 * ( Variance_wmRed + Variance_wmRed.' );
        [ ~, cholVariance_wmRed ] = NearestSPD( Variance_wmRed );

        Mean_wmRedMMedian_wmRed = Mean_wmRed - Median_wmRed;
        CholVariance_wmRed_Mean_wmRedMMedian_wmRed = cholVariance_wmRed * Mean_wmRedMMedian_wmRed;
        CholVariance_wmRed_Mean_wmRedMMedian_wmRed2 = CholVariance_wmRed_Mean_wmRedMMedian_wmRed.' * CholVariance_wmRed_Mean_wmRedMMedian_wmRed;

        if CholVariance_wmRed_Mean_wmRedMMedian_wmRed2 > eps && SkewLikelihood
            Zcheck_wmRed = ( Mean_wmRedMMedian_wmRed' * ano ) / realsqrt( CholVariance_wmRed_Mean_wmRedMMedian_wmRed2 );

            meanZcheck_wmRed = Zcheck_wmRed * CubatureWeights.';
            Zcheck_wmRed = Zcheck_wmRed - meanZcheck_wmRed;
            meanZcheck_wmRed2 = ( Zcheck_wmRed .* Zcheck_wmRed ) * CubatureWeights.';
            Zcheck_wmRed = Zcheck_wmRed / realsqrt( meanZcheck_wmRed2 );

            sZ3 = realpow( Zcheck_wmRed, 3 ) * CubatureWeights.';
            sZ4 = max( 3, realpow( Zcheck_wmRed, 4 ) * CubatureWeights.' );

            if DynamicNu
                tauno_nuno = LMFnlsq2( @( in ) CalibrateMomentsEST( in( 1 ), 4 + eps( 4 ) + exp( in( 2 ) ), Mean_wmRed, Median_wmRed, cholVariance_wmRed, sZ3, sZ4 ), [ 2; 0 ] );
                tauno = tauno_nuno( 1 );
                nuno = 4 + eps( 4 ) + exp( tauno_nuno( 2 ) );
            else
                tauno = LMFnlsq2( @( in ) CalibrateMomentsEST( in( 1 ), nuno, Mean_wmRed, Median_wmRed, cholVariance_wmRed, sZ3, [] ), 2 );
            end
        else
            tauno = Inf;

            if DynamicNu
                Zcheck_wmRed = ano;

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
            tauno = min( 100, tauno );
        end
        
        [ ~, wmno, deltaetano, CholPRRQno ] = CalibrateMomentsEST( tauno, nuno, Mean_wmRed, Median_wmRed, cholVariance_wmRed, [], [] );
        
        PersistentStateInternal_t = InvGetESTParametersFromVector( wmno, CholPRRQno, deltaetano, tauno, nuno, DynamicNu, SkewLikelihood );
        
    end
    
    p0 = PersistentStateInternal_t;
    f0 = ExpectedESTNLogPDF( p0, wmRed, CubatureWeights, Inf, DynamicNu, SkewLikelihood, nuno );
    
    if coder.target( 'MATLAB' )
        fminlbfgsOptions = struct( 'Display', 'iter', 'GradObj', 'on', 'GradConstr', true, 'GoalsExactAchieve', false, 'TolX', 1e-12, 'TolFun', 1e-12, 'MaxIter', Inf, 'MaxFunEvals', Inf, ...
            'HessUpdate', 'bfgs ', 'DiffMaxChange', 1e-1, 'DiffMinChange', 1e-8, 'OutputFcn', [], 'rho', 0.0100, 'sigma', 0.900, 'tau1', 3, 'tau2', 0.1, 'tau3', 0.5, 'StoreN', 20 );
    else
        fminlbfgsOptions = struct( 'Display', 'off', 'GradObj', 'on', 'GradConstr', true, 'GoalsExactAchieve', false, 'TolX', 1e-12, 'TolFun', 1e-12, 'MaxIter', Inf, 'MaxFunEvals', Inf, ...
            'HessUpdate', 'bfgs ', 'DiffMaxChange', 1e-1, 'DiffMinChange', 1e-8, 'OutputFcn', [], 'rho', 0.0100, 'sigma', 0.900, 'tau1', 3, 'tau2', 0.1, 'tau3', 0.5, 'StoreN', 20 );
    end
    [ pOpt, fOpt ] = fminlbfgs( @( p ) ExpectedESTNLogPDF( p, wmRed, CubatureWeights, f0 + 1, DynamicNu, SkewLikelihood, nuno ), p0, fminlbfgsOptions );
    if fOpt < f0
        PersistentStateInternal_t = pOpt;
    else
        PersistentStateInternal_t = p0;
    end
    
end
