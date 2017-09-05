function [ PersistentState, EndoSimulation, MeasurementSimulation ] = StochasticVolatilitySimulation( Parameters, PersistentState, InitialStates, ShockSequence, ~ )

    logit_mu = Parameters( 1 );
    phi = exp( -Parameters( 2 ) );
    phi = ( 1 - phi ) / ( 1 + phi );
    omega = exp( Parameters( 3 ) );
    rho = exp( -Parameters( 4 ) );
    rho = ( 1 - rho ) / ( 1 + rho );
    
    max_sigma = 1;
    
    if isempty( InitialStates )
        if isstruct( PersistentState ) && isfield( PersistentState, 'FinalEndo' ) && ( ~isempty( PersistentState.FinalEndo ) )
            CEndo = PersistentState.FinalEndo;
        else
            CEndo = [ logit_mu; 0 ];
        end
        T = size( ShockSequence, 2 );
        EndoSimulation = zeros( 2, T );
        
        for t = 1 : T
            CEndo = SimulationStep( CEndo, ShockSequence( :, t ), logit_mu, phi, omega, rho, max_sigma );
            EndoSimulation( :, t ) = CEndo;
        end
        
        PersistentState.FinalEndo = CEndo;
    else
        T = size( ShockSequence, 2 );
        EndoSimulation = zeros( 2, T );
        
        for t = 1 : T
            EndoSimulation( :, t ) = SimulationStep( InitialStates( :, t ), ShockSequence( :, t ), logit_mu, phi, omega, rho, max_sigma );
        end
    end
    
    if nargout > 2
        MeasurementSimulation = log( .40065867782212247344 + ( EndoSimulation( 2, : ) ./ max_sigma ).^ 2 ); % Maple to generate this number: f:=x->1/sqrt(2*Pi)*exp(-x^2/2); fsolve(int(log(c+x^2)*f(x),x=-infinity..infinity),c=0.5);
    end

end

function Endo = SimulationStep( Endo, Shocks, logit_mu, phi, omega, rho, max_sigma )
    L = [ 1, 0; rho, realsqrt( 1 - rho * rho ) ];
    CorrelatedShocks = L * Shocks;
    NewLogitSigma = ( 1 - phi ) * logit_mu + phi * Endo( 1 ) + omega * CorrelatedShocks( 1 );
    Endo = [ NewLogitSigma; ( max_sigma ./ ( 1 + exp( -NewLogitSigma ) ) ) * CorrelatedShocks( 2 ) ];
end
