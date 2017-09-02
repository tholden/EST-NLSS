% Copyright (c) 2013, John D'Errico & 2016, Tom Holden
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
% 

function [ Ahat, CholAhat ] = NearestSPD( A )
    % nearestSPD - the nearest (in Frobenius norm) Symmetric Positive Definite matrix to A
    % usage: Ahat = nearestSPD(A)
    %
    % From Higham: "The nearest symmetric positive semidefinite matrix in the
    % Frobenius norm to an arbitrary real matrix A is shown to be (B + H)/2,
    % where H is the symmetric polar factor of B=(A + A')/2."
    %
    % http://www.sciencedirect.com/science/article/pii/0024379588902236
    %
    % arguments: (input)
    %  A - square matrix, which will be converted to the nearest Symmetric
    %    Positive Definite Matrix.
    %
    % Arguments: (output)
    %  Ahat - The matrix chosen as the nearest SPD matrix to A.
    %  CholAhat - The real cholesky root of Ahat

    if isreal( A )
        [ Ahat, CholAhat ] = NearestSPDInternal( A );
    else
        Areal = real( A );
        [ ~, CholArealhat ] = NearestSPDInternal( Areal );
        Aimag = imag( A );
        Aimag = 0.5 * ( Aimag + Aimag.' );
        [ U, D ] = schur( Aimag, 'real' );
        d = diag( D );
        CholAhat = complex( CholArealhat );
        for j = 1 : length( d )
            if d( j ) ~= 0
                CholAhat = RealCholeskyUpdate( CholAhat, sqrt( d( j ) * 1i ) * U( :, j ), '+' );
            end
        end
        Ahat = CholAhat.' * CholAhat;
    end
        
end

function [ Ahat, CholAhat ] = NearestSPDInternal( A )

    ESTNLSSassert( nargin == 1, 'ESTNLSS:NearestSPD:Arguments', 'Exactly one argument must be provided to NearestSPD.' );
    ESTNLSSassert( all( isfinite( A(:) ) ), 'ESTNLSS:NearestSPD:NonFiniteInput', 'The input to NearestSPD must be finite.' );

    % test for a square matrix A
    [ r, c ] = size( A );
    ESTNLSSassert( r == c, 'ESTNLSS:NearestSPD:NonSquare', 'The input to NearestSPD must be a square matrix.' );
    
    if r == 1
        % A was scalar
        Ahat = max( real( A ), eps );
        CholAhat = realsqrt( Ahat );
        ESTNLSSassert( all( isfinite( Ahat(:) ) ), 'ESTNLSS:NearestSPD:NonFiniteOutputAhat', 'The Ahat output from NearestSPD was non-finite.' );
        ESTNLSSassert( all( isfinite( CholAhat(:) ) ), 'ESTNLSS:NearestSPD:NonFiniteOutputAhat', 'The cholAhat output from NearestSPD was non-finite.' );
        return
    end

    % symmetrize A into B
    Ahat = 0.5 * ( A + A' );

    [ CholAhat, p ] = chol( Ahat );

    if p ~= 0

        % Compute the symmetric polar factor of B. Call it H.
        % Clearly H is itself SPD.
        [ ~, Sigma, V ] = svd( Ahat );
        H = V * Sigma * V';

        % get Ahat in the above formula
        Ahat = 0.5 * ( Ahat + H );

        % ensure symmetry
        Ahat = 0.5 * ( Ahat + Ahat' );

        % test that Ahat is in fact PD. if it is not so, then tweak it just a bit.
        p = 1;
        for k = 0:100
            [ CholAhat, p ] = chol( Ahat );
            if p == 0
                break;
            end
            % Ahat failed the chol test. It must have been just a hair off,
            % due to floating point trash, so it is simplest now just to
            % tweak by adding a tiny multiple of an identity matrix.
            EigAhat =  eig( Ahat );
            IScale = abs( min( EigAhat ) );
            IScale = IScale + eps( IScale );
            IScale = max( IScale, eps( max( EigAhat ) ) );
            Ahat = Ahat + ( IScale * ( k .* k ) ) * eye( size( A ) );
        end
        ESTNLSSassert( p == 0, 'ESTNLSS:NearestSPD:Failure', 'Failed to find the nearest semi-positive definite matrix.' );

    end
    
    ESTNLSSassert( all( isfinite( Ahat(:) ) ), 'ESTNLSS:NearestSPD:NonFiniteOutputAhat', 'The Ahat output from NearestSPD was non-finite.' );
    ESTNLSSassert( all( isfinite( CholAhat(:) ) ), 'ESTNLSS:NearestSPD:NonFiniteOutputAhat', 'The cholAhat output from NearestSPD was non-finite.' );

end
