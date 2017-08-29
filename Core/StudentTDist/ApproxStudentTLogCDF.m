function log_y = ApproxStudentTLogCDF( x, nu )

    assert( numel( nu ) == 1, 'ESTNLSS:ApproxStudentTLogCDF:NuSize', 'ApproxStudentTLogCDF only supports univariate nu.' );
    assert( real( nu ) >= 0, 'ESTNLSS:ApproxStudentTLogCDF:NuSign', 'ApproxStudentTLogCDF requires nu to be weakly positive.' );
    assert( all( ~isnan( x(:) ) ), 'ESTNLSS:ApproxStudentTLogCDF:NaNInputX', 'ApproxStudentTLogCDF was passed a NaN input x.' );
    
    lognu = log( nu );
    if real( lognu ) <= -Inf
        log_y = log( 0.5 ) * ones( size( x ) );
        log_y( x <= -Inf ) = -Inf;
        log_y( x >= Inf ) = 0;
        return
    end
    
    log_y = zeros( size( x ) );
    
    SelInf = real( x ) >= Inf;
    
    SelNInf = real( x ) <= -Inf;
    log_y( SelNInf ) = -Inf;
    
    Done = SelInf | SelNInf;
    Remaining = find( ~Done );
    
    logH = log( 0.5 );
    
    xRemaining = x( Remaining );
    
    Pos_xRemaining = real( xRemaining ) > 0;
    xRemaining( Pos_xRemaining ) = -xRemaining( Pos_xRemaining );
    
    if real( nu ) < Inf
        
        SelSmallX = real( xRemaining .* xRemaining ) < real( nu );
        IdxSmallX = Remaining( SelSmallX );
        SelBigX = ~SelSmallX;
        IdxBigX = Remaining( SelBigX );
        
        Inu = 1 / nu;
        InuR2 = Inu * Inu;
        InuR3 = InuR2 * Inu;
        InuR4 = InuR2 * InuR2;
        InuR5 = InuR3 * InuR2;
        InuR6 = InuR3 * InuR3;
        
        cnu = [ +0;
                +0.3333333333333333333333 - 0.6666666666666666666667*Inu;
                +0.01111111111111111111111 + 0.3111111111111111111111*InuR2 -  0.1777777777777777777778*Inu;
                -0.0003527336860670194003527 -  0.2003527336860670194004*InuR3 + 0.1227513227513227513228*InuR2 -  0.01058201058201058201058*Inu;
                -0.00001763668430335097001764 +  0.1470194003527336860670*InuR4 - 0.09425044091710758377425*InuR3  + 0.009312169312169312169312*InuR2 + 0.0005643738977072310405644/ nu;
                +2.137779915557693335471e-6 -  0.1158163647052535941425*InuR5 + 0.07672064560953449842339*InuR4  - 0.008174870397092619314842*InuR3 - 0.0006840895729784618673508/ nu^2 + 0.00002992891881780770669660*Inu;
                +6.003533340394010588014e-9 + 0.09539614895523184764807*InuR6 -  0.06481074637688394302151*InuR5 +  0.007241589294499347409400*InuR4 +  0.0007510493295325746824865*InuR3 - 0.00003613709433815253921074/ nu^2 - 5.362460388915415370442e-6*Inu;
              ];
        
        t = xRemaining( SelSmallX );
        t = t(:).';
        
        tR2 = t .* t;
        tR4 = tR2 .* tR2;
        tR6 = tR4 .* tR2;
        tR8 = tR4 .* tR4;
        tR10 = tR6 .* tR4;
        tR12 = tR6 .* tR6;
        
        tR2OtR2Pnu = tR2 ./ ( tR2 + nu );
        
        cbetaln_nuO2_1O2 = cbetaln( nu * 0.5, 0.5 );
        
        log_y( IdxSmallX ) = log( 0.5 - exp( Shanks( bsxfun( @times, cnu, [ ones( size( t ) ); tR2; tR4; tR6; tR8; tR10; tR12 ] ) ) - cbetaln_nuO2_1O2 + 0.5 * log( tR2OtR2Pnu ) + ( 0.5 * nu - 1 ) * log1p( - tR2OtR2Pnu ) ) );
        
        InuP2 = 1 / ( 2 + nu );
        InuP2R2 = InuP2 * InuP2;
        InuP2R3 = InuP2R2 * InuP2;
        InuP2R4 = InuP2R2 * InuP2R2;
        InuP2R5 = InuP2R3 * InuP2R2;
        InuP2R6 = InuP2R3 * InuP2R3;
        InuP4 = 1 / ( 4 + nu );
        InuP4R2 = InuP4 * InuP4;
        InuP4R3 = InuP4R2 * InuP4;
        InuP6 = 1 / ( 6 + nu );
        InuP6R2 = InuP6 * InuP6;
        InuP8 = 1 / ( 8 + nu );
        InuP10 = 1 / ( 10 + nu );
        InuP12 = 1 / ( 12 + nu );
        
        Cnu = [ -lognu;
                -1+2*InuP2;
                +2.5-2*InuP2R2+8*InuP2-24*InuP4;
                -12.33333333333333333333+2.666666666666666666667*InuP2R3-16*InuP2R2+41*InuP2-288*InuP4+405*InuP6;
                +88.25-4*InuP2R4+32*InuP2R3-114*InuP2R2+239.5*InuP2-288*InuP4R2-2256*InuP4+9112.5*InuP6-8960*InuP8;
                -816.2+6.4*InuP2R5-64*InuP2R4+292*InuP2R3-807*InuP2R2+1516.166666666666666667*InuP2-6912*InuP4R2-10560*InuP4+123018.75*InuP6-334506.6666666666666667*InuP8+246093.75*InuP10;
                +9200.833333333333333333-10.66666666666666666667*InuP2R6+128*InuP2R5-712*InuP2R4+2440.666666666666666667*InuP2R3-5788.833333333333333333*InuP2R2+10131.82638888888888889*InuP2-4608*InuP4R3-95616*InuP4R2+17904*InuP4-82012.5*InuP6R2+1283495.625*InuP6-7470648.888888888888889*InuP8+13842773.4375*InuP10-8083152*InuP12;
              ];
        
        t = xRemaining( SelBigX );
        t = t(:).';
        
        tR2 = t .* t;
        ItR2 = 1 ./ tR2;
        ItR4 = ItR2 .* ItR2;
        ItR6 = ItR4 .* ItR2;
        ItR8 = ItR4 .* ItR4;
        ItR10 = ItR6 .* ItR4;
        ItR12 = ItR6 .* ItR6;
        
        nuOtR2Pnu = nu ./ ( tR2 + nu );
        
        log_y( IdxBigX ) = Shanks( bsxfun( @times, Cnu, [ ones( size( t ) ); ItR2; ItR4; ItR6; ItR8; ItR10; ItR12 ] ) ) - cbetaln_nuO2_1O2 + 0.5 * nu * log( nuOtR2Pnu ) - 0.5 * log1p( - nuOtR2Pnu );
        
        log_y( IdxBigX( real( ItR12 ) >= Inf ) ) = logH;
        
        log_y( IdxBigX( real( log_y( IdxBigX ) ) > real( logH ) ) ) = logH;
        
    else
        
        yRemaining = normcdf( xRemaining );
        SelSmallX = real( yRemaining ) > 0;
        IdxSmallX = Remaining( SelSmallX );
        SelBigX = ~SelSmallX;
        IdxBigX = Remaining( SelBigX );
        log_y( IdxSmallX ) = log( yRemaining( SelSmallX ) );
        
        t = xRemaining( SelBigX );
        t = t(:).';
        
        tR2 = t .* t;
        ItR2 = 1 ./ tR2;
        ItR4 = ItR2 .* ItR2;
        ItR6 = ItR4 .* ItR2;
        ItR8 = ItR4 .* ItR4;
        ItR10 = ItR6 .* ItR4;
        ItR12 = ItR6 .* ItR6;
        
        log_y( IdxBigX ) = Shanks( bsxfun( @times, [ 1; -1; 2.5; -12.33333333333333333333; 88.25; -816.2; 9200.833333333333333333 ], [ -.5*tR2-.9189385332046727417802-log(-t); ItR2; ItR4; ItR6; ItR8; ItR10; ItR12 ] ) );
        
        log_y( IdxBigX( real( ItR12 ) >= Inf ) ) = logH;
        
        log_y( IdxBigX( real( log_y( IdxBigX ) ) > real( logH ) ) ) = logH;
        
    end
    
    IdxPos = Remaining( Pos_xRemaining );
    log_y( IdxPos ) = log1p( -exp( log_y( IdxPos ) ) );
    
    assert( all( ~isnan( log_y(:) ) ), 'ESTNLSS:ApproxStudentTLogCDF:NaNOutputLogY', 'ApproxStudentTLogCDF returned a NaN output log_y.' );    
    
end
