function [f, B]=inbeta(z,a,b,N)
    %F=INBETA(X,P,Q) In complete beta function for complex x
    % there are branch cuts [-infinity,0] and [1 infinity] on the real axis
    % f = int(y.^(a-1).*(1-y).^(b-1),0..x);
    % B is the complete integral
    % The calculation is based in Muir's continued fraction expansion
    % Biometrika Vol. 22 284--297 (1930-31)
    % discussed in
    % L.A. Aroian Annals of Mathematical Statistics Vol. 12 218-223 (1941) 
    % and correection vol 30 1265 (1959)
    %
    % for large $z$ we use a Taylor expansion in powers of 1./z 
    %
    % The number of terms N can be adjusted
    %
    % Copyright Jim McElwaine 2010
    %
    % Call without arguments to run a test with random a and b in the
    % complex plane
    % The maximum error in the derivative is returned as a test
    
    if nargin<4
      N=[];
    end

    if isempty(N)
      N=20;
    end
    
    m0 = zeros( size( z ) );
    a = a + m0;
    b = b + m0;

    zmax = 5;
    s0 = (a+1)./(a+b+2);
    B = exp( cbetaln(a,b) );
    az = abs(z);
    f1 = find(real(z)<=s0 & az<zmax);
    f2 = find(real(z)>s0  & az<zmax);
    f3 = find(az>=zmax);

    f = complex(NaN(size(z)));
    f(f1) = inbeta3(z(f1),a(f1),b(f1),N);
    f(f2) = B(f2)-inbeta3(1-z(f2),b(f2),a(f2),N);
    f(f3) = inbeta2(z(f3),a(f3),b(f3),B(f3),N);

end

% Taylor series in z  
% expansion about 0
function f=inbeta1(z,a,b,N) %#ok<DEFNU>
    C=z.^a./a;
    f=C;
    for n=coder.unroll(1:N)
        C = z.*(C.*(n-b).*(n+a-1)./(n.*(n+a)));
        f = f+C;
        if max(abs(C))<tol
            break
        end
    end
end

% Taylor series in 1./z  
% expansion about infinity
function f=inbeta2(z,a,b,B,N)
    sz = sign(imag(z));
    C = -z.^(a+b-1)./(a+b-1).*exp(-1i.*b.*pi.*sz);
    z = 1./z;
    f = C;
    for n=coder.unroll(1:N)
        C = z.*C.*((n-b).*(n-a-b)./(n.*(n+1-a-b)));
        f = f+C;
    end
    f=f+B.*sin(pi.*b).*exp(1i.*pi.*a.*sz)./sin(pi.*a+pi.*b);
end

% Continued fraction expansion about z=0
function f=inbeta3(z,a,b,N)
    f=complex(zeros(size(z)));
    for k=coder.unroll(N:-1:1)
        f = k.*(b-k).*z./((a+2.*k-1).*(a+2.*k).*(1+f));
        j=k-1;
        f = -(a+j).*(a+b+j).*z./((a+2.*j).*(a+2.*j+1).*(1+f));
    end
    f=z.^a.*(1-z).^b./(a.*(1+f));
end

