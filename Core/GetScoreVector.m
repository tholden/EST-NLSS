function [ ObservationLikelihoods, PersistentState ] = GetScoreVector( InputParameters, Options, PersistentState )
    [ ~, PersistentState, ObservationLikelihoods ] = EstimationObjective( InputParameters, Options, PersistentState, false );
end
