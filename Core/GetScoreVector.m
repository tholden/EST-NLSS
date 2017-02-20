function [ ObservationLikelihoods, PersistentState ] = GetScoreVector( InputParameters, PersistentState, ObjectiveFunction )
    [ ~, PersistentState, ObservationLikelihoods ] = ObjectiveFunction( InputParameters, PersistentState );
end
