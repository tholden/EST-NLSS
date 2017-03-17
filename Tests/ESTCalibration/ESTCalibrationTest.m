clear all; %#ok<CLALL>
addpath ../../Core
addpath ../../Core/CholeskyUpdate/InbuiltImplementation
addpath ../../Core/StudentTDist
addpath ../../Core/ESTDist
addpath ../../Core/Utils

Order = 4;
N = 10;

xi = 10 * randn( N, 1 );
RootOmega = 0.1 * randn( N, N );
Omega = RootOmega * RootOmega';
[ Omega, cholOmega ] = NearestSPD( Omega );
delta = randn( N, 1 );
tau = 10 * randn;
nu = 4.5 + 4 * randn ^ 2; % Inf;

disp( 'tau, nu:' );
disp( [ tau, nu ] );

IntDim = N + 2;

log_tcdf_tau_nu = StudentTLogCDF( tau, nu );

if log_tcdf_tau_nu == 0
    IntDim = IntDim - 1;
    cholOmega = CholeskyUpdate( cholOmega, delta );
    Omega = cholOmega' * cholOmega;
    delta = zeros( N, 1 );
end

if isfinite( nu )
    [ Weights, pTmp, T ] = fwtpts( IntDim, Order );
    disp( 'T:' );
    disp( T );
    PhiN10 = normcdf( pTmp( end, : ) );
    if log_tcdf_tau_nu < 0
        N11Scaler = realsqrt( 0.5 * ( nu + 1 ) ./ gammaincinv( PhiN10, 0.5 * ( nu + 1 ), 'upper' ) );
    else
        N11Scaler = realsqrt( 0.5 * nu ./ gammaincinv( PhiN10, 0.5 * nu, 'upper' ) );
    end
end

if ~isfinite( nu ) || all( abs( N11Scaler - 1 ) <= realsqrt( eps ) )
    IntDim = IntDim - 1;
    [ Weights, pTmp, T ] = fwtpts( IntDim, Order );
    disp( 'T:' );
    disp( T );
else
    pTmp( end, : ) = [];
end

if log_tcdf_tau_nu < 0
    PhiN0 = normcdf( pTmp( end, : ) );
    pTmp( end, : ) = [];
    
    ICDFTmp = 1 - ( 1 - PhiN0 ) * exp( log_tcdf_tau_nu );
    FInvEST = zeros( size( ICDFTmp ) );
    FInvESTSelectLeft = ICDFTmp <= 0.5;
    FInvEST( FInvESTSelectLeft ) = StudentTInvLogCDF( reallog( ICDFTmp( FInvESTSelectLeft ) ), nu );
    FInvEST( ~FInvESTSelectLeft ) = -StudentTInvLogCDF( log_tcdf_tau_nu + reallog( 1 - PhiN0( ~FInvESTSelectLeft ) ), nu );
    
    tpdfRatio = exp( StudentTLogPDF( tau, nu ) - log_tcdf_tau_nu );
    
    ICDFTmp = 1 - 0.5 * exp( log_tcdf_tau_nu );
    if ICDFTmp <= 0.5
        MedT = StudentTInvLogCDF( reallog( ICDFTmp ), nu );
    else
        MedT = -StudentTInvLogCDF( log_tcdf_tau_nu - 0.693147180559945, nu ); % log( 0.5 ) = -0.693147180559945
    end
    
    N11Scaler = N11Scaler .* realsqrt( ( nu + FInvEST .^ 2 ) / ( 1 + nu ) );
else
    FInvEST = zeros( size( Weights ) );
    tpdfRatio = 0;
    MedT = 0;
end

ESTPoints = bsxfun( @plus, cholOmega' * bsxfun( @times, pTmp, N11Scaler ) + bsxfun( @times, delta, FInvEST ), xi );

mu = sum( bsxfun( @times, ESTPoints, Weights ), 2 );
DemeanedESTPoints = bsxfun( @minus, ESTPoints, mu );
Weighted_DemeanedESTPoints = bsxfun( @times, DemeanedESTPoints, Weights );

Sigma = DemeanedESTPoints * Weighted_DemeanedESTPoints';
[ Sigma, cholSigma ] = NearestSPD( Sigma );

lambda = ESTPoints( :, 1 );

tauTtau = tau * tau;
OPtauTtauDnu = 1 + tauTtau / nu;
if isfinite( nu )
    nuOnuM1 = nu / ( nu - 1 );
else
    nuOnuM1 = 1;
end
if log_tcdf_tau_nu < 0
    ET1 = nuOnuM1 * OPtauTtauDnu * tpdfRatio;
else
    ET1 = 0;
end
xiAlt = mu - delta * ET1;

lambdaAlt = xi + delta * MedT;

disp( 'xi, xiAlt:' );
disp( [ xi, xiAlt ] );

disp( 'lambda, lambdaAlt:' );
disp( [ lambda, lambdaAlt ] );

cholSigma_muMlambda = cholSigma * ( mu - lambda );

Zcheck = ( ( mu - lambda )' * DemeanedESTPoints ) / realsqrt( cholSigma_muMlambda' * cholSigma_muMlambda );

meanZcheck = Zcheck * Weights';
meanZcheck2 = Zcheck.^2 * Weights';

disp( 'EZ, EZ^2:' );
disp( [ meanZcheck, meanZcheck2 ] );

Zcheck = Zcheck - meanZcheck;
meanZcheck2 = Zcheck.^2 * Weights';
Zcheck = Zcheck / realsqrt( meanZcheck2 );

[ fZcheck, xiZcheck ] = ksdensity( Zcheck, linspace( min( Zcheck ), max( Zcheck ), 2000 ), 'NumPoints', 2000, 'Weights', Weights );
fZcheck = max( 0, fZcheck );
fZcheck = fZcheck / sum( fZcheck );
idx1Zcheck = find( fZcheck / max( fZcheck ) > 0.005, 1 );
idx2Zcheck = find( fZcheck / max( fZcheck ) > 0.005, 1, 'last' );
plot( xiZcheck( idx1Zcheck : idx2Zcheck ), fZcheck( idx1Zcheck : idx2Zcheck ) );

sZ3 = Zcheck.^3 * Weights';
sZ4 = max( 3, Zcheck.^4 * Weights' );

disp( 'EZ^3, EZ^3:' );
disp( [ sZ3, sZ4 ] );

[ resid, xiHat, deltaHat, cholOmegaHat ] = CalibrateMomentsEST( tau, nu, mu, lambda, cholSigma, sZ3, sZ4 );

disp( 'at truth:' );
disp( 'resid:' );
disp( resid' );
disp( 'xi comparison:' );
disp( [ xi, xiHat ] );
disp( 'delta comparison:' );
disp( [ delta, deltaHat ] );
disp( 'diag( cholOmega ) comparison:' );
disp( [ diag( cholOmega ), diag( cholOmegaHat ) ] );

Estim4 = LMFnlsq2( @( in ) CalibrateMomentsEST( in( 1 ), 4 + eps( 4 ) + exp( in( 2 ) ), mu, lambda, cholSigma, sZ3, sZ4 ), [ min( 10, tau ); reallog( min( 100, nu ) - 4 ) ] );
Estim3 = LMFnlsq2( @( in ) CalibrateMomentsEST( in( 1 ), nu, mu, lambda, cholSigma, sZ3, [] ), min( 10, tau ) );

Estim4( 2 ) = 4 + eps( 4 ) + exp( Estim4( 2 ) );

disp( 'Estim4 Estim3 Truth:' );
disp( [ Estim4( 1 ), Estim3, tau; Estim4( 2 ), nu, nu ] );

[ resid, xiHat, deltaHat, cholOmegaHat ] = CalibrateMomentsEST( Estim4( 1 ), Estim4( 2 ), mu, lambda, cholSigma, sZ3, sZ4 );

disp( 'at Estim4:' );
disp( 'resid:' );
disp( resid' );
disp( 'xi comparison:' );
disp( [ xi, xiHat ] );
disp( 'delta comparison:' );
disp( [ delta, deltaHat ] );
disp( 'diag( cholOmega ) comparison:' );
disp( [ diag( cholOmega ), diag( cholOmegaHat ) ] );

[ resid, xiHat, deltaHat, cholOmegaHat ] = CalibrateMomentsEST( Estim3( 1 ), nu, mu, lambda, cholSigma, sZ3, sZ4 );

disp( 'at Estim3:' );
disp( 'resid:' );
disp( resid' );
disp( 'xi comparison:' );
disp( [ xi, xiHat ] );
disp( 'delta comparison:' );
disp( [ delta, deltaHat ] );
disp( 'diag( cholOmega ) comparison:' );
disp( [ diag( cholOmega ), diag( cholOmegaHat ) ] );
