function [ Parameters, PersistentState ] = RunSmoothing( Parameters, Options, PersistentState )

    CorePath = [ fileparts( which( 'RunSmoothing' ) ) '/Core/' ];
    addpath( CorePath );

    NumParameters = size( Parameters, 1 );
    NumObservables = size( Options.Data, 1 );
    
    Options = SetDefaultOptions( Options );
    
    InputParameters = [ Parameters; bsxfun( @plus, log( Options.InitialMEStd ), zeros( NumObservables, 1 ) ) ];
        
    EstimatedNu = ~Options.NoTLikelihood && ~Options.DynamicNu;
    if EstimatedNu
        InputParameters( end + 1 ) = log( Options.InitialNu );
    end    
    
    rmpath( CorePath );

end
