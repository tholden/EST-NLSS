clear all; %#ok<CLALL>
addpath ../../Core
addpath ../../Core/CholeskyUpdate/InbuiltImplementation
addpath ../../Core/StudentTDist
addpath ../../Core/ESTDist
addpath ../../Core/Utils

N = 1;
FilterCubatureDegree = 20;
AllowTailEvaluations = true;

xi = 0; % 10 * randn( N, 1 );
RootOmega = 10; % 0.1 * randn( N, N );
Omega = RootOmega * RootOmega';
[ Omega, cholOmega ] = NearestSPD( Omega );
delta = 1; % randn( N, 1 );
tau = 1; % 10 * randn;
nu = 10; % 4.5 + 4 * randn ^ 2;

disp( 'tau, nu:' );
disp( [ tau, nu ] );

[ Weights, ESTPoints, NCubaturePoints, ET1, MedT ] = GetESTCubaturePoints( xi, Omega, delta, tau, nu, FilterCubatureDegree, eps ^ 0.375, AllowTailEvaluations );

disp( 'ET1, MedT:' );
disp( [ ET1, MedT ] );

mu = sum( bsxfun( @times, ESTPoints, Weights ), 2 );
muAlt = xi + delta * ET1;

lambda = ESTPoints( :, 1 );
lambdaAlt = xi + delta * MedT;

disp( 'mu, muAlt:' );
disp( [ mu, muAlt ] );

disp( 'lambda, lambdaAlt:' );
disp( [ lambda, lambdaAlt ] );

DemeanedESTPoints = bsxfun( @minus, ESTPoints, mu );
Weighted_DemeanedESTPoints = bsxfun( @times, DemeanedESTPoints, Weights );

Sigma = DemeanedESTPoints * Weighted_DemeanedESTPoints';
[ Sigma, cholSigma ] = NearestSPD( Sigma );

cholSigma_muMlambda = cholSigma * ( mu - lambda );

Zcheck = ( ( mu - lambda )' * DemeanedESTPoints ) / realsqrt( cholSigma_muMlambda' * cholSigma_muMlambda );

meanZcheck = Zcheck * Weights';
meanZcheck2 = ( Zcheck .* Zcheck ) * Weights';

disp( 'EZ, EZ^2:' );
disp( [ meanZcheck, meanZcheck2 ] );

Zcheck = Zcheck - meanZcheck;
meanZcheck2 = ( Zcheck .* Zcheck ) * Weights';
Zcheck = Zcheck / realsqrt( meanZcheck2 );

[ fZcheck, xiZcheck ] = ksdensity( Zcheck, linspace( min( Zcheck ), max( Zcheck ), 2000 ), 'NumPoints', 2000, 'Weights', Weights );
fZcheck = max( 0, fZcheck );
fZcheck = fZcheck / sum( fZcheck );
idx1Zcheck = find( fZcheck / max( fZcheck ) > 0.005, 1 );
idx2Zcheck = find( fZcheck / max( fZcheck ) > 0.005, 1, 'last' );
plot( xiZcheck( idx1Zcheck : idx2Zcheck ), fZcheck( idx1Zcheck : idx2Zcheck ) );

sZ3 = realpow( Zcheck, 3 ) * Weights';
sZ4 = max( 3, realpow( Zcheck, 4 ) * Weights' );

disp( 'EZ^3, EZ^4:' );
disp( [ sZ3, sZ4 ] );

Estim4 = LMFnlsq2( @( in ) CalibrateMomentsEST( in( 1 ), 4 + eps( 4 ) + exp( in( 2 ) ), mu, lambda, cholSigma, sZ3, sZ4 ), [ min( 10, tau ); reallog( min( 100, nu ) - 4 ) ] );
Estim3 = LMFnlsq2( @( in ) CalibrateMomentsEST( in( 1 ), nu, mu, lambda, cholSigma, sZ3, [] ), min( 10, tau ) );

Estim4( 2 ) = 4 + eps( 4 ) + exp( Estim4( 2 ) );

disp( 'Estim4 Estim3 Truth:' );
disp( [ Estim4( 1 ), Estim3, tau; Estim4( 2 ), nu, nu ] );

fprintf( '\n' );

disp( 'at truth:' );
DisplayResults( tau, nu, mu, lambda, cholSigma, sZ3, sZ4, xi, delta, cholOmega );

disp( 'at Estim4:' );
DisplayResults( Estim4( 1 ), Estim4( 2 ), mu, lambda, cholSigma, sZ3, sZ4, xi, delta, cholOmega );

disp( 'at Estim3:' );
DisplayResults( Estim3( 1 ), nu, mu, lambda, cholSigma, sZ3, sZ4, xi, delta, cholOmega );
