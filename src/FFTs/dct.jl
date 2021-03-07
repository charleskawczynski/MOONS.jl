"""
    dct2!(a, fft_type::FFTType=Radix2Type())

Discrete cosine transform of type II

Y = DCT2(X) returns the discrete cosine transform of X,
based on the staggered-grid definition
           N
   Y(k) =  ∑  X(j) cos (π*(k-1)*(j-1/2)/N)
          j=1
The vector Y is the same size as X and contains the
discrete cosine transform coefficients.

Y = DCT2(X,N) pads or truncates the vector X to length N
before transforming.

If X is a matrix, the DCT2 operation is applied to each
column. This transform can be inverted using IDCT2.

Verified correct within machine precision for 2D array input 2/28/2021

Traditional dct2:
          N-1
   Y(k) =  ∑  X(j) cos (π*(j+1/2)*k/N)
          j=0

"""
function dct2!(x, fft_type=Radix2Type())
    x .= dct2(x, fft_type)
end
function dct2(a, fft_type::FFTType=Radix2Type())
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
      aa[1:size(a,1),:] = a;
    else
      aa = a[1:n,:];
    end
    y = zeros(2*n, m);
    y[1:n,:] .= aa;
    i = Complex(0,1);
    e = 0.5 * exp.(-i*FT(0.5) * π * (0:n-1) ./ n);
    y[n+1:2*n,:] .= reverse(aa[1:n,:]; dims=1);
    yy = fft(y, 1);
    b = zeros(n,m);
    for l = 1:m
        b[1:n,l] .= real.(yy[1:n,l] .* e);
    end
    if isreal(a)
        b = real.(b);
    end
    if do_trans
        b = b'
    end
    return b
end
