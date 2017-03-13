function [ y, log_y ] = StudentTCDF( x, nu )

    assert( numel( nu ) == 1, 'ESTNLSS:StudentTCDF:NuSize', 'StudentTCDF only supports univariate nu.' );
    assert( nu > 0, 'ESTNLSS:StudentTCDF:NuSign', 'StudentTCDF requires nu to be strictly positive.' );
    assert( all( ~isnan( x(:) ) ), 'ESTNLSS:StudentTCDF:NaNInputX', 'StudentTCDF was passed a NaN input x.' );
    
    y = zeros( size( x ) );
    log_y = zeros( size( x ) );
    
    SelInf = x == Inf;
    y( SelInf ) = 1;
    
    SelNInf = x == -Inf;
    y( SelNInf ) = 0;
    log_y( SelNInf ) = -Inf;
    
    SelNaN = isnan( x );
    y( SelNaN ) = NaN;
    log_y( SelNaN ) = NaN;
    
    Done = SelInf | SelNInf | SelNaN;
    Remaining = find( ~Done );
    
    if nu < Inf && nu * reallog( nu ) < Inf
        y( Remaining ) = tcdf( x( Remaining ), nu );
        SelGood = Remaining( y( Remaining ) > 0 );
        SelBad = Remaining( y( Remaining ) == 0 );
        log_y( SelGood ) = reallog( y( SelGood ) );
        t = x( SelBad );
        kappa = 1 / nu;
        tkappa = t * kappa;
        tkappa2 = tkappa .* tkappa;
        tkappaM2 = 1 ./ tkappa2;
        tkappaM4 = tkappaM2 .* tkappaM2;
        tkappaM6 = tkappaM4 .* tkappaM2;
        tkappaM8 = tkappaM4 .* tkappaM4;
        log_y( SelBad ) = -reallog(nu)+0.5*nu*reallog(nu)-nu*reallog(-t)-betaln(0.5, 0.5*nu) + ( (-kappa-1)/(4*kappa+2) ) * tkappaM2 + ( kappa*(6*kappa^3+12*kappa^2+7*kappa+1)/(4*(2*kappa+1)^2*(4*kappa+1)) ) * tkappaM4 + ( -kappa^2*(60*kappa^5+170*kappa^4+176*kappa^3+80*kappa^2+15*kappa+1)/(6*(2*kappa+1)^3*(4*kappa+1)*(6*kappa+1)) ) * tkappaM6 + ( (3360*kappa^11+12880*kappa^10+20252*kappa^9+16850*kappa^8+8000*kappa^7+2200*kappa^6+346*kappa^5+29*kappa^4+kappa^3)/(8*(2*kappa+1)^4*(4*kappa+1)^2*(6*kappa+1)*(8*kappa+1)) ) * tkappaM8;
        y( SelBad ) = exp( log_y( SelBad ) );
        
        Remaining2 = SelBad( ~isfinite( log_y( SelBad ) ) | ~isfinite( y( SelBad ) ) );
    else
        Remaining2 = Remaining;
    end

    y( Remaining2 ) = normcdf( x( Remaining2 ) );
    SelGood = Remaining2( y( Remaining2 ) > 0 );
    SelBad = Remaining2( y( Remaining2 ) == 0 );
    log_y( SelGood ) = reallog( y( SelGood ) );
    t = x( SelBad );
    t2 = t .* t;
    tM2 = 1 ./ t2;
    tM4 = tM2 .* tM2;
    tM6 = tM4 .* tM2;
    tM8 = tM4 .* tM4;
    log_y( SelBad ) = -.500000000000000000*t2-.918938533204672742-reallog(-t)-tM2+2.50000000000000000*tM4-12.3333333333333333*tM6+88.2500000000000000*tM8;
    y( SelBad ) = exp( log_y( SelBad ) );
    
    assert( all( isfinite( y(:) ) ), 'ESTNLSS:StudentTCDF:NonFiniteOutputY', 'StudentTCDF returned a non-finite output y.' );
    assert( all( ~isnan( log_y(:) ) ), 'ESTNLSS:StudentTCDF:NaNOutputLogY', 'StudentTCDF returned a NaN output log_y.' );    
    
end
