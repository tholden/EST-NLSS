function DisplayResults( tau, nu, mu, lambda, cholSigma, sZ3, sZ4, xi, delta, cholOmega )

    [ resid, xiHat, deltaHat, cholOmegaHat, Z3, Z4 ] = CalibrateMomentsEST( tau, nu, mu, lambda, cholSigma, sZ3, sZ4 );

    fprintf( '\n' );
    disp( 'EZ^3 EZ^4:' );
    disp( [ Z3, Z4 ] );
    disp( 'resid:' );
    disp( resid' );
    disp( 'xi comparison:' );
    disp( [ xi, xiHat ] );
    disp( 'delta comparison:' );
    disp( [ delta, deltaHat ] );
    disp( 'diag( cholOmega ) comparison:' );
    disp( [ diag( cholOmega ), diag( cholOmegaHat ) ] );
    fprintf( '\n' );

end
