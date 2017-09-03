function Success = ESTNLSSSetup

    if verLessThan( 'matlab', '9.1' )
        disp( 'EST-NLSS compilation requires MATLAB r2016b or later.' );
        return
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
    
    disp( 'Building KalmanStepInternal. This may take some time.' );
    
    try
        KalmanStepInternal_script;
    catch Error
        disp( 'Error building KalmanStepInternal:' );
        DisplayError( Error );
        SuccessInternal = false;
    end
    
    disp( 'Building EstimationObjectiveInternal. This may take some time.' );
    
    try
        EstimationObjectiveInternal_script;
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
