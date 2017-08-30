function [ log_y, Dlog_y ] = ExpectedESTLogPDF( p, X, W )

    n = size( X, 1 );
    
    [ xi, CholOmega, delta, tau, nu ] = GetESTParametersFromVector( p, n );
    
    
end
