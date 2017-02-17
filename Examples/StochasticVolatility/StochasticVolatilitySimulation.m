function [ PersistentState, EndoSimulation, MeasurementSimulation ] = StochasticVolatilitySimulation( Parameters, PersistentState, InitialStates, ShockSequence, ~ )

    mu = Parameters( 1 );
    phi = Parameters( 2 );
    omega = Parameters( 3 );
    rho = Parameters( 4 );
    
    if isempty( InitialStates )
    else
    end

end

function NewEndo = SimulationStep( OldEndo, Shocks, mu, phi, omega, rho )
    L = [ 1, 0; rho, sqrt( 1 - rho * rho ) ];
    CorrelatedShocks = L * Shocks;
    NewLogSigma = ( 1 - phi ) * mu + phi * OldEndo( 1 ) + omega * CorrelatedShocks( 1 );
    NewEndo = [ NewLogSigma; exp( NewLogSigma ) * CorrelatedShocks( 2 ) ];
end
