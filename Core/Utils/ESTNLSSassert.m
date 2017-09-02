function ESTNLSSassert( cond, msgID, msg, varargin )
    if ~cond
        ESTNLSSerror( msgID, msg, varargin{:} );
    end
end
