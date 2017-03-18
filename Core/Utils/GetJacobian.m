function Jacobian = GetJacobian( f, x, nf )
    nx = length( x );
    Jacobian = NaN( nf, nx );
    if ~all( isfinite( f( x ) ) )
        return;
    end
    sreps = sqrt( eps );
    srsreps = sqrt( sreps );
    parfor i = 1 : nx
        WarningState = warning( 'off', 'all' );
        xi = x( i );
        h = max( [ sreps, abs( srsreps * xi ) ] );
        while true
            if h < eps
                h = eps;
            end
            try
                Jacobian( :, i ) = ( f( SetElement( x, i, xi + h ) ) - f( SetElement( x, i, xi - h ) ) ) / ( 2 * h ); %#ok<PFBNS>
            catch
            end
            if all( isfinite( Jacobian( :, i ) ) ) || h == eps
                break;
            else
                h = 0.5 * h;
            end
        end
        warning( WarningState );
    end
end

function x = SetElement( x, i, xi )
    x( i ) = xi;
end
