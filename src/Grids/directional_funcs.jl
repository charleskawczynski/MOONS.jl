##### directional funcs

export orthog, orthogs

orthogs(::Val{1}) = (2,3)
orthogs(::Val{2}) = (3,1)
orthogs(::Val{3}) = (1,2)
orthogs(i::Int) = orthogs(Val(i))

orthog(::Val{3}, ::Val{2}) = 1
orthog(::Val{2}, ::Val{3}) = 1
orthog(::Val{3}, ::Val{1}) = 2
orthog(::Val{1}, ::Val{3}) = 2
orthog(::Val{1}, ::Val{2}) = 3
orthog(::Val{2}, ::Val{1}) = 3
orthog(i::Int, j::Int) = orthog(Val(i), Val(j))
