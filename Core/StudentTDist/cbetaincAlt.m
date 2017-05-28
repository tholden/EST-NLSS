% // <copyright file="Beta.cs" company="Math.NET">
% // Math.NET Numerics, part of the Math.NET Project
% // http://numerics.mathdotnet.com
% // http://github.com/mathnet/mathnet-numerics
% //
% // Copyright (c) 2009-2014 Math.NET
% //
% // Permission is hereby granted, free of charge, to any person
% // obtaining a copy of this software and associated documentation
% // files (the "Software"), to deal in the Software without
% // restriction, including without limitation the rights to use,
% // copy, modify, merge, publish, distribute, sublicense, and/or sell
% // copies of the Software, and to permit persons to whom the
% // Software is furnished to do so, subject to the following
% // conditions:
% //
% // The above copyright notice and this permission notice shall be
% // included in all copies or substantial portions of the Software.
% //
% // THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% // EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
% // OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% // NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
% // HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
% // WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% // FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
% // OTHER DEALINGS IN THE SOFTWARE.
% // </copyright>
% 
% // <contribution>
% //    Cephes Math Library, Stephen L. Moshier
% //    ALGLIB 2.0.1, Sergey Bochkanov
% // </contribution>

function y = cbetaincAlt( x, a, b )

    y = zeros( size( x ) );

    a = a + y;
    b = b + y;
    
    bt = y;
    
    bt( ( x > 0 ) & ( x < 1 ) ) = exp( cgammaln( a + b ) - cgammaln( a ) - cgammaln( b ) + ( a .* log( x ) ) + ( b .* log( 1 - x ) ) );
    
	ST = x >= ( a + 1 ) ./ ( a + b + 2 );

	% Continued fraction representation
	eeps = eps;
	fpmin = eps( 0 ) ./ eeps;

    x( ST ) = 1 - x( ST );
    swap = a;
    a( ST ) = b( ST );
    b( ST ) = swap( ST );

	qab = a + b;
	qap = a + 1;
	qam = a - 1;
	c = 1;
	d = 1 - (qab.*x./qap);

	d = cRectify( d, fpmin );

	d = 1./d;
	h = d;

    coder.unroll( );
	for m = 1 : 140
	
        m2 = 2 .* m;
        
		aa = m.*(b - m).*x./((qam + m2).*(a + m2));
		d = 1 + (aa.*d);

		d = cRectify( d, fpmin );

		c = 1 + (aa./c);
		c = cRectify( c, fpmin );

		d = 1./d;
		h = h .* d.*c;
		aa = -(a + m).*(qab + m).*x./((a + m2).*(qap + m2));
		d = 1 + (aa.*d);

		d = cRectify( d, fpmin );

		c = 1 + (aa./c);

		c = cRectify( c, fpmin );

		d = 1./d;
		del = d.*c;
		h = h .* del;

		if all( abs( del - 1 ) <= eeps )
			break
		end
	end

    y( ST ) = 1 - (bt(ST).*h(ST)./a(ST));
    y( ~ST ) = bt(~ST).*h(~ST)./a(~ST);

end

function z = cRectify( z, fpmin )
    fz = isfinite( z );
    abs_z = abs( z );
    sel = abs_z < fpmin;
    z( sel ) = ( z( sel ) ./ abs_z( sel ) ) .* fpmin;
    z( ~isfinite( z ) & fz ) = fpmin;
end
