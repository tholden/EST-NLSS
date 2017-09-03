function out = LeftScalerOp( op, a, b )
    coder.inline( 'always' );
    if coder.target( 'MATLAB' )
        out = op( a, b );
    else
        out = coder.nullcopy( b );
        for i = 1 : numel( b )
            out( i ) = op( a, b( i ) );
        end
    end
end
