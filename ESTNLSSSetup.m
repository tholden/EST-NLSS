function Success = ESTNLSSSetup( DebugMode )

    if verLessThan( 'matlab', '9.1' )
        disp( 'EST-NLSS compilation requires MATLAB r2016b or later.' );
        return
    end
    
    if nargin < 1
        DebugMode = false;
    end

    CorePath = [ fileparts( which( 'ESTNLSSSetup' ) ) '/Core/' ];

    addpath( CorePath );
    addpath( [ CorePath 'Utils' ] );
    addpath( [ CorePath 'StudentTDist' ] );
    addpath( [ CorePath 'ESTDist' ] );

    if verLessThan( 'matlab', '9.2' )
        addpath( [ CorePath 'CholeskyUpdate/MImplementation/' ] );
    else
        addpath( [ CorePath 'CholeskyUpdate/InbuiltImplementation/' ] );
    end

    rehash;

    OldPath = cd( CorePath );
    
    SuccessInternal = true;

    disp( 'Building KalmanStepInternal. This may take several hours.' );

    try
        if DebugMode
            KalmanStepInternal_scriptDebug;
        else
            KalmanStepInternal_script;
        end
    catch Error
        disp( 'Error building KalmanStepInternal:' );
        DisplayError( Error );
        SuccessInternal = false;
    end

    disp( 'Building EstimationObjectiveInternal. This may take several hours.' );

    try
        if DebugMode
            EstimationObjectiveInternal_scriptDebug;
        else
            EstimationObjectiveInternal_script;
        end
    catch Error
        disp( 'Error building EstimationObjectiveInternal:' );
        DisplayError( Error );
        SuccessInternal = false;
    end

    cd( OldPath );
    
    if nargout > 0
        Success = SuccessInternal;
    end
    
end
