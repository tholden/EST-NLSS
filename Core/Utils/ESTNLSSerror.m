function ESTNLSSerror( msgID, msg, varargin )
    if coder.target( 'MATLAB' )
        warning( msgID, msg, varargin{:} );
        % error( msgID, msg, varargin{:} );
        % keyboard;
    else
        fprintf( 'Error: %s\n', msgID );
        fprintf( [ msg '\n' ], varargin{:} );
    end
end
