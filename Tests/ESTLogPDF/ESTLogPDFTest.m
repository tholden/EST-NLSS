clear all; %#ok<CLALL>
addpath ../../Core
addpath ../../Core/CholeskyUpdate/InbuiltImplementation
addpath ../../Core/StudentTDist
addpath ../../Core/ESTDist
addpath ../../Core/Utils

rng( 'default' );

N = 1;

xi = 10 * randn( N, 1 );
RootOmega = 0.1 * randn( N, N );
Omega = RootOmega * RootOmega';
[ Omega, CholOmega ] = NearestSPD( Omega );
delta = randn( N, 1 );
tau = 10 * randn;
nu = 2 + 4 * randn ^ 2;

disp( 'tau, nu:' );
disp( [ tau, nu ] );

disp( 'Integral of PDF using CholOmega:' );
disp( integral( @( x ) exp( ESTLogPDF( x, xi, CholOmega, delta, tau, nu, true ) ), -100, 100 ) );

disp( 'Integral of PDF using Omega:' );
disp( integral( @( x ) exp( ESTLogPDF( x, xi, Omega, delta, tau, nu, false ) ), -100, 100 ) );

N = 2;

xi = 10 * randn( N, 1 );
RootOmega = 0.1 * randn( N, N );
Omega = RootOmega * RootOmega';
[ Omega, CholOmega ] = NearestSPD( Omega );
delta = randn( N, 1 );
tau = 10 * randn;
nu = 2 + 4 * randn ^ 2;

disp( 'tau, nu:' );
disp( [ tau, nu ] );

disp( 'Integral of PDF using CholOmega:' );
disp( integral2( @( x1, x2 ) reshape( exp( ESTLogPDF( [ x1(:).'; x2(:).' ], xi, CholOmega, delta, tau, nu, true ) ), size( x1 ) ), -100, 100, -100, 100 ) );

disp( 'Integral of PDF using Omega:' );
disp( integral2( @( x1, x2 ) reshape( exp( ESTLogPDF( [ x1(:).'; x2(:).' ], xi, Omega, delta, tau, nu, false ) ), size( x1 ) ), -100, 100, -100, 100 ) );

N = 3;

xi = 10 * randn( N, 1 );
RootOmega = 0.1 * randn( N, N );
Omega = RootOmega * RootOmega';
[ Omega, CholOmega ] = NearestSPD( Omega );
delta = randn( N, 1 );
tau = 10 * randn;
nu = 2 + 4 * randn ^ 2;

disp( 'tau, nu:' );
disp( [ tau, nu ] );

disp( 'Integral of PDF using CholOmega:' );
disp( integral3( @( x1, x2, x3 ) reshape( exp( ESTLogPDF( [ x1(:).'; x2(:).'; x3(:).' ], xi, CholOmega, delta, tau, nu, true ) ), size( x1 ) ), -100, 100, -100, 100, -100, 100 ) );

disp( 'Integral of PDF using Omega:' );
disp( integral3( @( x1, x2, x3 ) reshape( exp( ESTLogPDF( [ x1(:).'; x2(:).'; x3(:).' ], xi, Omega, delta, tau, nu, false ) ), size( x1 ) ), -100, 100, -100, 100, -100, 100 ) );
