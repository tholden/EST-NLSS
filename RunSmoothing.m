function [ Parameters, PersistentState ] = RunSmoothing( Parameters, Options, PersistentState )

    CorePath = [ fileparts( which( 'RunSmoothing' ) ) '/Core/' ];
    addpath( CorePath );

    rmpath( CorePath );

end
