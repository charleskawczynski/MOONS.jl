"""
    idct2!(x, fft_type::FFTType = Radix2Type())

Inverse discrete cosine transform of type II.

X = IDCT2(Y) inverts the DCT2 transform, returning the
original vector if Y was obtained using Y = DCT2(X). It is based on
the staggered-grid definition
                          N
   X(j) = 2/N*( Y(1)/2 +  ∑ Y(k) cos (π*(k-1)*(j-1/2)/N) )
                         k=2

X = IDCT2(Y,N) pads or truncates the vector Y to length N
before transforming.

If Y is a matrix, the IDCT2 operation is applied to
each column.

                          N
   X(j) = 2/N*( Y(1)/2 +  ∑ Y(k) cos (π*(k-1)*(j-1/2)/N) )
                         k=2
"""
function idct2!(x, fft_type=Radix2Type())
    x .= idct2(x, fft_type)
end
function idct2(a, fft_type::FFTType=Radix2Type())
    FT = eltype(a)
    if min(size(a)...)==1
        if size(a,2)>1
            do_trans = true;
        else
            do_trans = false;
        end
        a = a[:];
    else
        do_trans = false;
    end
    n = size(a,1);
    m = size(a,2);
    # Pad or truncate a if necessary
    if size(a,1)<n
      aa = zeros(n,m);
      aa[1:size(a,1),:] .= a;
    else
      aa = a[1:n,:];
    end
    y = Array{Complex{FT},2}(undef, 2*n, m)
    fill!(y, 0)
    y[1:n,:] .= aa;
    i = Complex(0,1);
    e = 0.5 * exp.(-i*FT(0.5) * π * (0:n-1) ./ n);
    for l = 1:m
        y[1:n,l] .= y[1:n,l] .* e;
    end
    y[n+2:2*n,:] .= reverse(conj.(y[2:n,:]); dims=1);
    yy = fft(y, 1);
    b = (2/n) * yy[1:n,:];
    if isreal(a)
        b = real.(b)
    end
    if do_trans
        b = b'
    end
    return b
end
