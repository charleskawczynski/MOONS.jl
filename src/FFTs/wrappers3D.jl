using UnPack
using ..Fields: sweep
using ..Fields: pre_post_colons

#####
##### fft!
#####

dct3D!(f, dir::Int) =
    dct3D!(f, Val(dir), Radix2Type())

function dct3D!(f, ::Val{dim}, fft_type::FFTType=Radix2Type()) where {dim}
    @unpack Ipre, Ipost = pre_post_colons(f, dim)
    for i in 1:size(f, dim)
        fv = @view f[Ipre..., i, Ipost...]
        dct2!(fv, fft_type)
    end
end

idct3D!(f, dir::Int) =
    idct3D!(f, Val(dir), Radix2Type())

function idct3D!(f, ::Val{dim}, fft_type::FFTType=Radix2Type()) where {dim}
    @unpack Ipre, Ipost = pre_post_colons(f, dim)
    for i in 1:size(f, dim)
        fv = @view f[Ipre..., i, Ipost...]
        idct2!(fv, fft_type)
    end
end

