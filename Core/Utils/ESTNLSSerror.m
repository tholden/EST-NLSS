function ESTNLSSerror( msgID, msg, varargin )
    if coder.target( 'MATLAB' )
        warning( msgID, msg, varargin{:} );
    end
    % else
        fprintf( 'Warning: %s\n', msgID );
        fprintf( [ msg '\n' ], varargin{:} );
    % end
end
