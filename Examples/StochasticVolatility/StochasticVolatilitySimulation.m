function [ PersistentState, EndoSimulation, MeasurementSimulation ] = StochasticVolatilitySimulation( Parameters, PersistentState, InitialStates, ShockSequence, ~ )

    log_mu = Parameters( 1 );
    phi = exp( -Parameters( 2 ) );
    phi = ( 1 - phi ) / ( 1 + phi );
    omega = exp( Parameters( 3 ) );
    rho = exp( -Parameters( 4 ) );
    rho = ( 1 - rho ) / ( 1 + rho );
    
    if isempty( InitialStates )
        CEndo = [ log_mu; 0 ];
        T = size( ShockSequence, 2 );
        EndoSimulation = zeros( 2, T );
        
        for t = 1 : T
            CEndo = SimulationStep( CEndo, ShockSequence( :, t ), log_mu, phi, omega, rho );
            EndoSimulation( :, t ) = CEndo;
        end
    else
        T = size( ShockSequence, 2 );
        EndoSimulation = zeros( 2, T );
        
        for t = 1 : T
            EndoSimulation( :, t ) = SimulationStep( InitialStates( :, t ), ShockSequence( :, t ), log_mu, phi, omega, rho );
        end
    end
    
    if nargout > 2
        MeasurementSimulation = EndoSimulation( 2, : );
    end

end

function Endo = SimulationStep( Endo, Shocks, log_mu, phi, omega, rho )
    L = [ 1, 0; rho, sqrt( 1 - rho * rho ) ];
    CorrelatedShocks = L * Shocks;
    NewLogSigma = ( 1 - phi ) * log_mu + phi * Endo( 1 ) + omega * CorrelatedShocks( 1 );
    Endo = [ NewLogSigma; exp( NewLogSigma ) * CorrelatedShocks( 2 ) ];
end
